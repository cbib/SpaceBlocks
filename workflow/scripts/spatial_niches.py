"""
spatial_niches.py – Cross-sample spatial niche / domain identification (BANKSY)
================================================================================
Runs BANKSY jointly across all preprocessed samples to find spatial domains
("niches") that are common across the cohort, then writes one niche label per
cell back as a per-sample TSV that annotate_cells injects as ``spatial_niche``.

Flow (multi-sample BANKSY with Harmony, cell resolution):
  1. Concatenate all preprocessed per-sample h5ads.
  2. (Optional) subset to HVGs, then z-score (BANKSY convention).
  3. Stagger each sample's spatial coordinates side-by-side so the single
     BANKSY kNN graph never links spots across samples.
  4. initialize_banksy + generate_banksy_matrix  → neighbour-augmented matrix.
  5. PCA on the BANKSY matrix.
  6. Harmony across samples on the BANKSY PCA (spatially-informed integration).
  7. Leiden at niche_resolution → spatial_niche.
  8. (Optional) BANKSY spatial majority-vote refinement of the niche labels
     (refine), applied to freshly-computed OR reloaded/precomputed niches.

Outputs:
  spatial_niches/
  ├── spatial_niches_concatenated.h5ad   (light: lognorm X + niche + embeddings)
  ├── tsv/niche_{sample}.tsv             (barcode → spatial_niche, per sample)
  └── plots/                             (spatial maps, UMAPs, composition)
"""

import logging
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
from scipy import sparse

from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy_utils.umap_pca import pca_umap
from banksy_utils.refine_clusters import refine_once

# Shared composition-barplot helpers (scripts/ is on sys.path for script: rules)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composition_barplots import composition_pair, build_niche_palette  # noqa: E402


# ── Logging ──────────────────────────────────────────────────────────────────
log_handlers = [logging.StreamHandler(sys.stderr)]
if hasattr(snakemake, "log"):
    if snakemake.log.out:
        Path(snakemake.log.out).parent.mkdir(parents=True, exist_ok=True)
        log_handlers.append(logging.FileHandler(snakemake.log.out, mode="w"))
    if snakemake.log.err:
        Path(snakemake.log.err).parent.mkdir(parents=True, exist_ok=True)
        sys.stderr = open(snakemake.log.err, "w")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=log_handlers,
)
log = logging.getLogger("spatial_niches")


# ── Parameters ───────────────────────────────────────────────────────────────
adata_paths     = [str(p) for p in snakemake.input.adatas]
sample_ids      = list(snakemake.params.sample_ids)
SEED            = int(snakemake.params.random_seed)
LAMBDA          = float(snakemake.params.banksy_lambda)
NUM_NEIGHBOURS  = int(snakemake.params.num_neighbours)
MAX_M           = int(snakemake.params.max_m)
NBR_DECAY       = str(snakemake.params.nbr_weight_decay)
USE_HVG         = bool(snakemake.params.use_hvg)
N_TOP_GENES     = int(snakemake.params.n_top_genes)
N_PCS           = int(snakemake.params.n_pcs)
NICHE_RES       = float(snakemake.params.niche_resolution)
COMPUTE_UMAP    = bool(getattr(snakemake.params, "compute_umap", False))
CLUSTER_NN      = int(getattr(snakemake.params, "cluster_n_neighbours", 15))
N_ITERATIONS    = int(getattr(snakemake.params, "n_iterations", 2))
USE_PRECOMPUTED = bool(getattr(snakemake.params, "use_precomputed", False))
NICHE_DIR       = str(getattr(snakemake.params, "niche_dir", "") or "")
REFINE          = bool(getattr(snakemake.params, "refine", False))
REFINE_NEIGH    = int(getattr(snakemake.params, "refine_num_neigh", 6))
REFINE_ITERS    = int(getattr(snakemake.params, "refine_iterations", 1))
REFINE_AUTO     = bool(getattr(snakemake.params, "refine_auto", False))
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors
DPI             = int(getattr(snakemake.params, "dpi", 300))

out_concat   = str(snakemake.output.concatenated)
niche_tsvs   = [str(p) for p in snakemake.output.niche_tsvs]
plots_dir    = str(snakemake.output.plots_dir)

NICHE_KEY = "spatial_niche"
COORD_KEYS = ("x_pixel", "y_pixel", "coord_xy")   # BANKSY reads obsm[coord_keys[2]]


def harmony_embedding(emb, obs, batch_key, seed):
    """Harmony-correct a PCA embedding directly (not scanpy's wrapper, which
    hard-codes a Z_corr transpose that breaks across harmonypy versions).
    Returns an (n_obs, n_pcs) array; on any degenerate output the uncorrected
    embedding is returned so the rule never crashes."""
    emb = np.asarray(emb, dtype=np.float32)
    n, k = emb.shape
    try:
        ho = harmonypy.run_harmony(emb, obs, [batch_key], random_state=seed)
        Z = np.asarray(ho.Z_corr)
        if Z.ndim == 2 and Z.shape == (k, n):          # old convention (d, N)
            Z = Z.T
        elif Z.ndim == 2 and Z.shape == (n, k):        # new convention (N, d)
            pass
        else:
            raise ValueError(f"unexpected Harmony Z_corr shape {Z.shape}")
        return np.ascontiguousarray(Z, dtype=np.float32)
    except Exception as e:
        log.warning("Harmony unavailable/degenerate (%s); using uncorrected PCA.", e)
        return emb


def stagger_coordinates(adata, sample_key, coord_keys):
    """Stagger per-sample spatial coordinates the way the BANKSY multi-sample
    notebook does: left-align each sample to x-min 0, then offset each sample
    along x by sample_index × (1.5 × the largest per-sample x-extent), leaving
    y unchanged. The staggered coords are written BOTH to obs[coord_keys[0]] /
    obs[coord_keys[1]] (used for plotting) AND to obsm[coord_keys[2]] — the
    latter is what initialize_banksy actually reads to build the spatial graph,
    so a single graph never links cells across samples. obsm['spatial'] is left
    untouched for per-sample plotting in native coordinates."""
    x_key, y_key, xy_key = coord_keys
    sp = np.asarray(adata.obsm["spatial"], dtype=float)
    df = pd.DataFrame({"x": sp[:, 0], "y": sp[:, 1],
                       "sample": np.asarray(adata.obs[sample_key].values)})
    df["x"] = df.groupby("sample")["x"].transform(lambda v: v - v.min())
    global_max_x = float(df["x"].max()) * 1.5
    # preserve the concat sample order for left-to-right layout
    order = list(pd.unique(df["sample"]))
    sample_no = pd.Categorical(df["sample"], categories=order).codes
    df["x"] = df["x"] + sample_no * global_max_x
    adata.obs[x_key] = df["x"].values
    adata.obs[y_key] = df["y"].values
    adata.obsm[xy_key] = np.vstack([df["x"].values, df["y"].values]).T
    log.info("  staggered %d samples along x (global_max_x=%.1f)",
             len(order), global_max_x)


def _spatial_scatter(adata, sample_key, color_key, out_path, dpi, color_map):
    """Per-sample spatial domain maps (scatter of obsm['spatial'], no image
    dependency), one panel per sample, using the shared niche colour map."""
    samples = list(pd.unique(adata.obs[sample_key]))
    cats = list(pd.Categorical(adata.obs[color_key]).categories)
    colour = {c: color_map.get(str(c), "#cccccc") for c in cats}
    n = len(samples)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow),
                             squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for i, s in enumerate(samples):
        ax = axes[i // ncol][i % ncol]
        ax.set_axis_on()
        m = (adata.obs[sample_key] == s).values
        sp = np.asarray(adata.obsm["spatial"])[m]
        vals = adata.obs[color_key].values[m]
        for c in cats:
            mm = vals == c
            if mm.any():
                ax.scatter(sp[mm, 0], sp[mm, 1], s=2, linewidths=0,
                           color=colour[c], rasterized=True)
        ax.set_title(str(s), fontsize=9)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=5,
                          color=colour[c], label=str(c)) for c in cats]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, title="Spatial niche", fontsize=7, title_fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _resolve_precomputed():
    """Return {sample: tsv_path} with a niche TSV for EVERY sample (external
    niche_dir first, else the rule's own output dir). SAFEGUARD: if any sample's
    TSV is missing, ERROR OUT — when use_precomputed is requested, silently
    recomputing would defeat the reproducibility intent, so we fail loudly rather
    than fall back to computing."""
    src, missing = {}, []
    for sid in sample_ids:
        cands = []
        if NICHE_DIR:
            cands.append(os.path.join(NICHE_DIR, f"niche_{sid}.tsv"))
        cands.append(os.path.join(tsv_dir, f"niche_{sid}.tsv"))  # pipeline output
        found = next((c for c in cands if os.path.isfile(c)), None)
        if found is None:
            missing.append(sid)
        else:
            src[sid] = found
    if missing:
        raise FileNotFoundError(
            f"use_precomputed=True but niche TSVs are MISSING for "
            f"{len(missing)}/{len(sample_ids)} sample(s): {', '.join(missing)} "
            f"(searched {NICHE_DIR or tsv_dir}). Aborting — set use_precomputed: "
            f"false to compute, or provide a complete set of niche_{{sample}}.tsv.")
    log.info("Precomputed niche TSVs resolved for all %d samples.", len(sample_ids))
    return src


try:
    log.info("=" * 70)
    log.info("Spatial niche identification (BANKSY) across %d samples", len(adata_paths))
    log.info("  lambda=%.2f  num_neighbours=%d  max_m=%d  decay=%s",
             LAMBDA, NUM_NEIGHBOURS, MAX_M, NBR_DECAY)
    log.info("  use_hvg=%s  n_pcs=%d  niche_resolution=%.2f  seed=%d",
             USE_HVG, N_PCS, NICHE_RES, SEED)
    log.info("=" * 70)

    os.makedirs(plots_dir, exist_ok=True)
    tsv_dir = os.path.dirname(niche_tsvs[0])
    os.makedirs(tsv_dir, exist_ok=True)

    # ── 1. Load + concatenate ────────────────────────────────────────────
    adatas = []
    for path, sid in zip(adata_paths, sample_ids):
        log.info("  Loading %s …", sid)
        a = sc.read_h5ad(path)
        if "sample" not in a.obs.columns:
            a.obs["sample"] = sid
        if "cell_id" not in a.obs.columns:
            a.obs["cell_id"] = a.obs_names.astype(str)
        adatas.append(a)

    # preprocess_umap already sets obs["sample"] = sample_id per cell, so we do
    # NOT pass label= (which would clash with that column). keys + index_unique
    # only disambiguate obs_names across samples.
    adata = sc.concat(adatas, join="inner", keys=sample_ids, index_unique="-")
    del adatas
    log.info("Concatenated: %d cells, %d genes", adata.n_obs, adata.n_vars)

    if "spatial" not in adata.obsm:
        raise ValueError("obsm['spatial'] missing after concatenation.")

    # ── 2. Stagger coordinates (both modes: BANKSY graph in compute mode,
    #       spatial plots in both) ───────────────────────────────────────
    log.info("Staggering per-sample spatial coordinates …")
    stagger_coordinates(adata, "sample", COORD_KEYS)

    # ── 3. Niche labels: reload precomputed, or compute with BANKSY ──────
    if USE_PRECOMPUTED:
        niche_source = _resolve_precomputed()   # ERRORS if any TSV is missing
        log.info("use_precomputed=True → reloading niche labels, skipping BANKSY/Leiden")
        labels = pd.Series(index=adata.obs_names, dtype=object)
        uncovered = {}
        for sid, src in niche_source.items():
            ndf = pd.read_csv(src, sep="\t", index_col=0)
            col = NICHE_KEY if NICHE_KEY in ndf.columns else (
                ndf.columns[0] if len(ndf.columns) else None)
            m = (adata.obs["sample"] == sid).values
            ids = adata.obs["cell_id"].values[m]
            if col is None or ndf.shape[0] == 0:
                uncovered[sid] = (int(m.sum()), int(m.sum()))      # empty TSV
                continue
            vals = ndf[col].reindex(ids)
            n_miss = int(vals.isna().sum())
            if n_miss:
                uncovered[sid] = (n_miss, int(m.sum()))
            labels.loc[adata.obs_names[m]] = vals.values
            log.info("  %s ← %s", sid, src)
        # SAFEGUARD: every cell of every sample must be covered, else ERROR OUT.
        if uncovered:
            detail = "; ".join(f"{s}: {n}/{tot} cells uncovered"
                               for s, (n, tot) in uncovered.items())
            raise ValueError(
                "use_precomputed=True but the niche TSVs do not fully cover all "
                f"cells — {detail}. Aborting (the TSV may be empty or out of sync "
                "with the current cells). Recompute with use_precomputed: false, or "
                "supply complete niche TSVs.")
        labels = labels.astype(str)
        try:
            raw_cats = sorted(labels.unique().tolist(), key=lambda x: int(x))
        except ValueError:
            raw_cats = sorted(labels.unique().tolist())
        adata.obs[NICHE_KEY] = pd.Categorical(
            labels.values, categories=[str(c) for c in raw_cats], ordered=True)
        log.info("Reloaded spatial_niche → %d niches", adata.obs[NICHE_KEY].nunique())
    else:
        # NB: adata.X is log-normalized (preprocess_umap: normalize_total +
        # log1p; raw counts in layers["raw_counts"]). The default "seurat"
        # flavor expects log-normalized data, so it is correct here.
        bank_input = adata.copy()
        if USE_HVG:
            log.info("Selecting %d HVGs (flavor=seurat, on log-normalized X) …",
                     N_TOP_GENES)
            sc.pp.highly_variable_genes(bank_input, n_top_genes=N_TOP_GENES,
                                        flavor="seurat")
            bank_input = bank_input[:, bank_input.var["highly_variable"]].copy()
            log.info("  Using %d HVGs", bank_input.n_vars)
        else:
            log.info("use_hvg=False → using all %d genes", bank_input.n_vars)

        # Memory guard: BANKSY matrix is DENSE and (2 + max_m)·n_genes wide.
        feats = bank_input.n_vars * (2 + MAX_M)
        est_gb = bank_input.n_obs * feats * 4 / 1e9
        log.info("BANKSY matrix ~%d cells × %d features (~%.1f GB dense, float32)",
                 bank_input.n_obs, feats, est_gb)
        if est_gb > 80:
            log.warning("  Projected BANKSY matrix is large (~%.0f GB). Consider "
                        "use_hvg=true / smaller n_top_genes if this OOMs.", est_gb)

        # ── BANKSY neighbour-augmented matrix (uses staggered coord_xy) ──
        log.info("Initialising BANKSY + building neighbour-augmented matrix …")
        banksy_dict = initialize_banksy(
            bank_input, coord_keys=COORD_KEYS,
            num_neighbours=NUM_NEIGHBOURS, nbr_weight_decay=NBR_DECAY, max_m=MAX_M,
            plt_edge_hist=False, plt_nbr_weights=False, plt_agf_angles=False,
            plt_theta=False,
        )
        banksy_dict, _ = generate_banksy_matrix(
            bank_input, banksy_dict, lambda_list=[LAMBDA], max_m=MAX_M, verbose=False,
        )

        # ── PCA (official banksy_utils.pca_umap) ─────────────────────────
        log.info("PCA (%d PCs) on the BANKSY matrix …", N_PCS)
        pca_umap(banksy_dict, pca_dims=[N_PCS], plt_remaining_var=False, add_umap=False)
        banksy_adata = banksy_dict[NBR_DECAY][LAMBDA]["adata"]
        adata.obsm["X_pca_banksy"] = np.asarray(
            banksy_adata.obsm[f"reduced_pc_{N_PCS}"], dtype=np.float32)

        # ── Harmony across samples on the BANKSY PCA ─────────────────────
        log.info("Harmony integration across samples …")
        adata.obsm["X_pca_harmony_banksy"] = harmony_embedding(
            banksy_adata.obsm[f"reduced_pc_{N_PCS}"], adata.obs, "sample", SEED)
        del bank_input, banksy_dict, banksy_adata

        # ── Clustering: scanpy approximate-NN graph + igraph Leiden ──────
        log.info("Building kNN graph (n_neighbors=%d) on the Harmony embedding …",
                 CLUSTER_NN)
        sc.pp.neighbors(adata, use_rep="X_pca_harmony_banksy",
                        n_neighbors=CLUSTER_NN, random_state=SEED)

        def _leiden(res, key):
            sc.tl.leiden(adata, resolution=res, key_added=key, random_state=SEED,
                         flavor="igraph", n_iterations=N_ITERATIONS, directed=False)

        log.info("Leiden (igraph, %d iterations) @ resolution %.2f …",
                 N_ITERATIONS, NICHE_RES)
        _leiden(NICHE_RES, NICHE_KEY)
        raw_cats = sorted(adata.obs[NICHE_KEY].unique().tolist(), key=lambda x: int(x))
        adata.obs[NICHE_KEY] = pd.Categorical(
            adata.obs[NICHE_KEY].astype(str),
            categories=[str(c) for c in raw_cats], ordered=True)
        log.info("spatial_niche @ resolution %.2f → %d niches",
                 NICHE_RES, adata.obs[NICHE_KEY].nunique())

        # Integrated UMAP (opt-in; viz-only, slow on ~1M cells)
        if COMPUTE_UMAP:
            try:
                import umap
                log.info("Computing integrated UMAP on the Harmony embedding …")
                adata.obsm["X_umap_banksy"] = umap.UMAP(
                    random_state=SEED).fit_transform(adata.obsm["X_pca_harmony_banksy"])
            except Exception as e:
                log.warning("UMAP skipped (%s)", e)
        else:
            log.info("compute_umap=False → skipping integrated UMAP")

    # ── 4. Optional BANKSY spatial refinement ────────────────────────────
    # Spatial majority-vote smoothing (BANKSY refine_once): each cell adopts its
    # spatial-neighbourhood majority niche when >50% of neighbours agree. Reads
    # the STAGGERED coord_xy, so neighbours are within-sample (per-sample
    # smoothing). Applies to BOTH freshly-computed and reloaded/precomputed
    # niches. Raw labels are kept in obs["{NICHE_KEY}_raw"].
    if REFINE:
        log.info("Spatial refinement (BANKSY majority vote, num_neigh=%d, %s) …",
                 REFINE_NEIGH, "auto" if REFINE_AUTO else f"{REFINE_ITERS} iter")
        adata.obs[f"{NICHE_KEY}_raw"] = adata.obs[NICHE_KEY].copy()
        cur = list(adata.obs[NICHE_KEY].astype(str))
        max_iter = 50 if REFINE_AUTO else max(1, REFINE_ITERS)
        for it in range(max_iter):
            new, _, _ = refine_once(adata, cur, None, COORD_KEYS,
                                    num_neigh=REFINE_NEIGH)
            changed = int(np.sum(np.asarray(cur) != np.asarray(new)))
            cur = list(new)
            ratio = changed / max(1, len(cur))
            log.info("  refine iter %d: %d cells changed (%.2f%%)",
                     it + 1, changed, 100 * ratio)
            if REFINE_AUTO and ratio < 0.005:
                break
            if not REFINE_AUTO and (it + 1) >= REFINE_ITERS:
                break
        try:
            ref_cats = sorted(set(cur), key=lambda x: int(x))
        except ValueError:
            ref_cats = sorted(set(cur))
        adata.obs[NICHE_KEY] = pd.Categorical(
            cur, categories=[str(c) for c in ref_cats], ordered=True)
        n0 = adata.obs[f"{NICHE_KEY}_raw"].nunique()
        n1 = adata.obs[NICHE_KEY].nunique()
        log.info("Refinement complete: %d → %d niches "
                 "(raw kept in obs['%s_raw'])", n0, n1, NICHE_KEY)

    # ── Per-sample niche TSVs (barcode → spatial_niche) ──────────────────
    log.info("Writing per-sample niche TSVs …")
    tsv_by_sample = {os.path.basename(p): p for p in niche_tsvs}
    for sid in sample_ids:
        m = (adata.obs["sample"] == sid).values
        df = pd.DataFrame(
            {NICHE_KEY: adata.obs[NICHE_KEY].values[m].astype(str)},
            index=adata.obs["cell_id"].values[m],
        )
        df.index.name = "barcode"
        out_p = tsv_by_sample.get(f"niche_{sid}.tsv",
                                  os.path.join(tsv_dir, f"niche_{sid}.tsv"))
        df.to_csv(out_p, sep="\t")
        log.info("  %s: %d cells, %d niches", sid, df.shape[0],
                 df[NICHE_KEY].nunique())

    # ── Plots ────────────────────────────────────────────────────────────
    log.info("Plots …")
    # ONE niche → colour mapping reused across every niche plot (per-sample
    # spatial maps, combined staggered map, UMAPs, composition barplots) so a
    # niche is the same colour everywhere.
    niche_cfg = (ANNOTATION_COLORS.get(NICHE_KEY, {})
                 if isinstance(ANNOTATION_COLORS, dict) else {})
    cats = list(adata.obs[NICHE_KEY].cat.categories)
    niche_palette = build_niche_palette(cats, niche_cfg)
    adata.uns[f"{NICHE_KEY}_colors"] = [niche_palette[str(c)] for c in cats]

    try:
        _spatial_scatter(adata, "sample", NICHE_KEY,
                         os.path.join(plots_dir, "spatial_niches_per_sample.png"),
                         DPI, niche_palette)
    except Exception as e:
        log.warning("  spatial scatter failed: %s", e)

    # Staggered-coordinate layout coloured by sample (verifies staggering)
    try:
        xy = adata.obsm[COORD_KEYS[2]]
        codes = pd.Categorical(adata.obs["sample"]).codes
        fig, ax = plt.subplots(figsize=(min(40, 4 * len(sample_ids)), 5))
        ax.scatter(xy[:, 0], xy[:, 1], c=codes, cmap="tab20", s=1,
                   linewidths=0, rasterized=True)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_axis_off()
        ax.set_title("Staggered coordinates (coloured by sample)", fontsize=10)
        fig.savefig(os.path.join(plots_dir, "staggered_coordinates_by_sample.png"),
                    dpi=DPI, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        log.warning("  staggered-coords plot failed: %s", e)

    # Combined staggered domain map coloured by niche (all samples, one view)
    try:
        xy = adata.obsm[COORD_KEYS[2]]
        colors = adata.obs[NICHE_KEY].astype(str).map(niche_palette).values
        fig, ax = plt.subplots(figsize=(min(40, 4 * len(sample_ids)), 5))
        ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=1,
                   linewidths=0, rasterized=True)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_axis_off()
        ax.set_title("BANKSY spatial niches (staggered, all samples)", fontsize=10)
        handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=5,
                              color=niche_palette[str(c)], label=str(c)) for c in cats]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
                  frameon=False, title="Niche", fontsize=6, title_fontsize=7,
                  ncol=max(1, len(cats) // 20 + 1))
        fig.savefig(os.path.join(plots_dir, "spatial_niches_staggered.png"),
                    dpi=DPI, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        log.warning("  combined staggered domain plot failed: %s", e)

    if "X_umap_banksy" in adata.obsm:
        for color, fname in [(NICHE_KEY, "umap_by_niche.png"),
                             ("sample", "umap_by_sample.png")]:
            try:
                sc.pl.embedding(adata, basis="X_umap_banksy", color=color, size=3,
                                frameon=False, show=False,
                                title=f"BANKSY (Harmony) – {color}")
                plt.savefig(os.path.join(plots_dir, fname), dpi=DPI,
                            bbox_inches="tight")
                plt.close()
            except Exception as e:
                log.warning("  UMAP (%s) failed: %s", color, e)

    # niche composition (by sample, and by region if present)
    try:
        composition_pair(adata, "sample", NICHE_KEY, niche_palette, plots_dir,
                         "niche_by_sample", group_label="Sample",
                         cat_label="Niche", dpi=DPI)
        has_regions = ("region_annotation" in adata.obs.columns
                       and adata.obs["region_annotation"].nunique() > 1
                       and not all(adata.obs["region_annotation"] == "Unlabeled"))
        if has_regions:
            composition_pair(adata, "region_annotation", NICHE_KEY, niche_palette,
                             plots_dir, "niche_by_region", group_label="Region",
                             cat_label="Niche", dpi=DPI)
    except Exception as e:
        log.warning("  niche composition barplots failed: %s", e)

    # ── Save light concatenated object ───────────────────────────────────
    log.info("Saving concatenated → %s", out_concat)
    Path(out_concat).parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_concat)

    log.info("Spatial niche identification complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
