"""
qc_sweep.py — OPTIONAL pre-filtering QC diagnostic (modular core, module 1).
============================================================================
Reads the VALIDATED, unfiltered, segmented h5ad (core contract input) and makes
the consequences of candidate QC thresholds explicit, so the final analysis.*
thresholds can be chosen deliberately before run_upstream. It does NOT filter,
NOT cluster and does NOT write an h5ad.

Per-feature thresholds are LISTS (e.g. min_genes=[100, 200]) resolved per sample
upstream (shared `default` + optional `per_sample` overrides). For every cutoff
it shows WHERE on the tissue the affected cells land.

Outputs (per sample)
--------------------
  * qc_violins.png            : counts / genes / %mito violins (grouped by region
                                if present, strip on, scatter rasterised) with
                                every candidate threshold drawn as a line.
  * qc_remove_low_genes.png   : per-parameter "landing" figure — leftmost panel =
  * qc_remove_low_counts.png    whole dataset, then one panel per threshold showing
  * qc_remove_high_counts.png   the removed cells on the tissue. Two rows when an
  * qc_remove_high_mito.png     ingest reference is given (top = by cell type,
                                bottom = plain), one row otherwise.
  * qc_features_spatial.png   : continuous spatial composite of the three metrics.
  * qc_joint_scatter.png      : n_genes vs total_counts coloured by %mito, with
                                thresholds drawn (doublet / debris view).
  * qc_thresholds_summary.tsv : per-feature ranges + per-threshold n removed /
                                retained and the percentile each threshold lands at.
  * qc_ingest_removed.tsv     : (if ingest) cell type × threshold → removed; else stub.

Mito is computed only when genes match the prefix(es); spatial plots use the
embedded image when present and otherwise scatter on coordinates.
"""
import logging
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import numpy as np
import pandas as pd
import scanpy as sc

# ── Logging ──────────────────────────────────────────────────────────────────
log_handlers = [logging.StreamHandler(sys.stderr)]
if hasattr(snakemake, "log"):
    if snakemake.log.out:
        Path(snakemake.log.out).parent.mkdir(parents=True, exist_ok=True)
        log_handlers.append(logging.FileHandler(snakemake.log.out, mode="w"))
    if snakemake.log.err:
        Path(snakemake.log.err).parent.mkdir(parents=True, exist_ok=True)
        sys.stderr = open(snakemake.log.err, "w")
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S", handlers=log_handlers)
log = logging.getLogger("qc_sweep")


try:
    h5ad_path  = str(snakemake.input.h5ad)
    sample_id  = str(snakemake.params.sample_id)
    thresholds = dict(snakemake.params.thresholds or {})
    mito_prefix = tuple(snakemake.params.mito_prefix)
    DPI = int(snakemake.params.dpi)
    ingest_enabled = bool(getattr(snakemake.params, "ingest_enabled", False))
    ref_label_key  = str(getattr(snakemake.params, "ref_label_key", "cell_type"))

    out = snakemake.output
    out_dir = os.path.dirname(str(out.violins_png))
    os.makedirs(out_dir, exist_ok=True)

    th_genes = thresholds.get("min_genes") or []
    th_clow  = thresholds.get("min_counts") or []
    th_chigh = thresholds.get("max_counts") or []
    th_mt    = thresholds.get("max_pct_mt") or []

    log.info("=" * 70)
    log.info("QC sweep (diagnostic, no filtering): sample=%s", sample_id)
    log.info("  min_genes=%s  min_counts=%s  max_counts=%s  max_pct_mt=%s",
             th_genes, th_clow, th_chigh, th_mt)
    log.info("=" * 70)

    adata = sc.read_h5ad(h5ad_path)
    adata.obs["sample"] = sample_id
    log.info("  loaded %d cells × %d genes (unfiltered)", adata.n_obs, adata.n_vars)

    # QC metrics (robust mito)
    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    mito_available = bool(adata.var["mt"].any())
    sc.pp.calculate_qc_metrics(adata, qc_vars=(["mt"] if mito_available else []),
                               percent_top=None, log1p=False, inplace=True)
    if not mito_available:
        th_mt = []
        log.warning("  no mito genes (prefix %s) — %%mito panels skipped", mito_prefix)

    region_key = ("region_annotation"
                  if ("region_annotation" in adata.obs.columns
                      and adata.obs["region_annotation"].nunique() > 1) else None)

    # One spatial conditional, reused everywhere. Always pass an explicit spot_size
    # (cell centroids are far smaller than the scalefactor-derived Visium spot size,
    # so without it the dots render invisibly and only the tissue image shows); add
    # library_id to draw the dots over the embedded image when one is present.
    spatial_kw = {"spot_size": 20}
    if isinstance(adata.uns.get("spatial"), dict) and adata.uns["spatial"]:
        spatial_kw["library_id"] = list(adata.uns["spatial"].keys())[0]
    has_coords = "spatial" in adata.obsm

    # Optional ingest overlay: transfer cell-type labels from the SAME reference the
    # downstream ingest_ref rule uses (config: ingest_ref + ingest_ref_label_key), so
    # the "what am I removing" panels can be broken down by cell type. Runs the same
    # sc.tl.ingest recipe as ingest_ref; the contract X is raw counts, so a copy of the
    # query is normalised/log1p'd on shared genes before projection. Diagnostic only —
    # never written back to the contract.
    ingest_vals, ingest_palette = None, {}
    ingest_ref_path = getattr(snakemake.input, "ingest_ref", None)
    if isinstance(ingest_ref_path, (list, tuple)):
        ingest_ref_path = ingest_ref_path[0] if ingest_ref_path else None
    if ingest_enabled and ingest_ref_path and os.path.exists(str(ingest_ref_path)):
        try:
            log.info("  ingest overlay: transferring '%s' from %s",
                     ref_label_key, ingest_ref_path)
            adata_ref = sc.read_h5ad(str(ingest_ref_path))
            if ref_label_key not in adata_ref.obs.columns:
                raise ValueError(f"label key '{ref_label_key}' not in reference obs")
            adata_ref.obs[ref_label_key] = adata_ref.obs[ref_label_key].astype("category")
            shared = adata.var_names.intersection(adata_ref.var_names)
            if len(shared) < 100:
                raise ValueError(f"only {len(shared)} shared genes — too few for ingest")
            # normalise a query copy (contract X = raw counts) into the reference space
            q = adata[:, shared].copy()
            sc.pp.normalize_total(q, target_sum=1e4)
            sc.pp.log1p(q)
            r = adata_ref[:, shared].copy()
            sc.pp.pca(r)
            sc.pp.neighbors(r)
            sc.tl.umap(r)
            sc.tl.ingest(q, r, obs=ref_label_key)
            ingest_vals = pd.Series(q.obs[ref_label_key].astype(str).values,
                                    index=adata.obs_names)
            base = (list(matplotlib.colormaps["tab20"].colors)
                    + list(matplotlib.colormaps["tab20b"].colors))
            ingest_palette = {c: matplotlib.colors.to_hex(base[i % len(base)])
                              for i, c in enumerate(sorted(ingest_vals.unique()))}
            log.info("  ingest labels transferred (%d types, %d shared genes)",
                     len(ingest_palette), len(shared))
            del q, r, adata_ref
        except Exception as e:
            log.warning("  ingest overlay failed (skipping): %s", e)
            ingest_vals = None

    tc = adata.obs["total_counts"].values
    ng = adata.obs["n_genes_by_counts"].values
    mt = adata.obs["pct_counts_mt"].values if mito_available else None

    # ── 1. Violins — composite of two rows: top WITH per-cell observations
    #        (stripplot, rasterised), bottom WITHOUT (clean distributions). Too many
    #        points obscure the threshold lines, so both views are shown together.
    #        Candidate thresholds are drawn on both rows.
    vpanels = [("total_counts", "Total counts", th_clow + th_chigh),
               ("n_genes_by_counts", "Genes per cell", th_genes)]
    if mito_available:
        vpanels.append(("pct_counts_mt", "% mitochondrial", th_mt))
    ncols = len(vpanels)
    fig, axes = plt.subplots(2, ncols, figsize=(6.0 * ncols, 10), squeeze=False)
    for r, (strip, row_tag) in enumerate([(True, "with observations"),
                                          (False, "distribution only")]):
        for c, (col, lab, lines) in enumerate(vpanels):
            ax = axes[r][c]
            try:
                sc.pl.violin(adata, col, groupby=region_key, ax=ax, show=False,
                             stripplot=strip, jitter=0.4, size=1, rotation=45)
                if strip:                              # rasterise the per-cell strip only
                    for art in ax.collections:
                        if isinstance(art, mcoll.PathCollection):
                            art.set_rasterized(True)
            except Exception as e:
                log.warning("  violin %s (%s) failed: %s", col, row_tag, e)
            for t in lines:
                ax.axhline(t, ls="--", lw=0.8, color="#e41a1c")
            ax.set_title(f"{lab}\n({row_tag})")
            ax.set_xlabel("")
    fig.suptitle(f"Unfiltered QC — {sample_id}  "
                 f"({adata.n_obs:,} cells, {adata.n_vars:,} genes)"
                 + ("" if mito_available else "  [no MT- genes]"),
                 fontsize=13, fontweight="bold", y=1.005)
    fig.tight_layout()
    fig.savefig(str(out.violins_png), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    log.info("  saved %s", out.violins_png)

    # ── 2. Per-parameter landing figures ─────────────────────────────────
    def make_landing(label, vals, op, ths, out_path):
        """One file per parameter: leftmost column = whole dataset, then one
        column per threshold; top row by ingest type (if available), bottom row
        plain (removed in red, kept grey)."""
        if not ths or not has_coords:
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.text(0.5, 0.5,
                    "no spatial coordinates" if not has_coords else f"{label}: no thresholds",
                    ha="center", va="center")
            ax.set_axis_off()
            fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
            plt.close(fig)
            return

        nrows = 2 if ingest_vals is not None else 1
        ncols = 1 + len(ths)
        fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 5.0 * nrows),
                                 squeeze=False)

        def draw(ax, color, title, **extra):
            try:
                sc.pl.spatial(adata, color=color, ax=ax, show=False, title=title,
                              **spatial_kw, **extra)
            except Exception as e:
                log.warning("  landing '%s' failed: %s", title, e)
                ax.set_title(f"{title} (failed)")
                ax.axis("off")

        r = 0
        if ingest_vals is not None:
            adata.obs["_qc"] = pd.Categorical(ingest_vals.astype(str))
            adata.uns["_qc_colors"] = [ingest_palette[c]
                                       for c in adata.obs["_qc"].cat.categories]
            draw(axes[r][0], "_qc", "all cells (ingest)",
                 legend_fontsize=5, na_in_legend=False)
            for j, t in enumerate(ths):
                mask = (vals < t) if op == "below" else (vals > t)
                title = f"{label} {t}\n{int(mask.sum()):,} removed"
                if mask.sum() == 0:
                    draw(axes[r][j + 1], None, title)
                    continue
                tmp = pd.Series(np.where(mask, ingest_vals.astype(str), None),
                                index=adata.obs_names)
                cats = [c for c in ingest_palette if c in set(tmp.dropna())]
                adata.obs["_qc"] = pd.Categorical(tmp, categories=cats)
                adata.uns["_qc_colors"] = [ingest_palette[c] for c in cats]
                draw(axes[r][j + 1], "_qc", title, na_color="#e8e8e8",
                     na_in_legend=False, legend_fontsize=5)
            r += 1

        draw(axes[r][0], None, "all cells")
        for j, t in enumerate(ths):
            mask = (vals < t) if op == "below" else (vals > t)
            title = f"{label} {t}\n{int(mask.sum()):,} removed"
            if mask.sum() == 0:
                draw(axes[r][j + 1], None, title)
                continue
            adata.obs["_qc"] = pd.Categorical(
                pd.Series(np.where(mask, "removed", None), index=adata.obs_names),
                categories=["removed"])
            adata.uns["_qc_colors"] = ["#e41a1c"]
            draw(axes[r][j + 1], "_qc", title, na_color="#e8e8e8", na_in_legend=False)

        if "_qc" in adata.obs:
            del adata.obs["_qc"]
        adata.uns.pop("_qc_colors", None)
        fig.suptitle(f"{label} — {sample_id}", fontsize=13, fontweight="bold", y=1.01)
        fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        log.info("  saved %s", out_path)

    make_landing("n_genes <", ng, "below", th_genes, str(out.low_genes_png))
    make_landing("counts <", tc, "below", th_clow, str(out.low_counts_png))
    make_landing("counts >", tc, "above", th_chigh, str(out.high_counts_png))
    make_landing("%mito >", mt, "above", th_mt, str(out.high_mito_png))

    # ── 3. Continuous QC features in space ───────────────────────────────
    feats = [("n_genes_by_counts", "Genes per cell"), ("total_counts", "Total counts")]
    if mito_available:
        feats.append(("pct_counts_mt", "% mitochondrial"))
    if has_coords:
        fig, axes = plt.subplots(1, len(feats), figsize=(7 * len(feats), 6),
                                  gridspec_kw={"wspace": 0.3})
        if len(feats) == 1:
            axes = [axes]
        for ax, (col, lab) in zip(axes, feats):
            try:
                sc.pl.spatial(adata, color=col, ax=ax, show=False, title=lab,
                              cmap="viridis", **spatial_kw)
            except Exception as e:
                log.warning("  feature spatial %s failed: %s", col, e)
                ax.set_title(f"{lab} (failed)")
                ax.axis("off")
        fig.suptitle(f"QC features in space — {sample_id}",
                     fontsize=13, fontweight="bold", y=1.02)
    else:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "no spatial coordinates", ha="center", va="center")
        ax.set_axis_off()
    fig.savefig(str(out.features_png), dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    # ── 4. Joint n_genes × total_counts coloured by %mito (+ thresholds) ─
    fig, ax = plt.subplots(figsize=(7, 6))
    sca = ax.scatter(ng, tc, c=(mt if mito_available else "#377eb8"), s=2, alpha=0.4,
                     cmap=("viridis" if mito_available else None),
                     linewidths=0, rasterized=True)
    if mito_available:
        fig.colorbar(sca, ax=ax).set_label("% mitochondrial")
    for t in th_genes:
        ax.axvline(t, ls="--", lw=0.8, color="#e41a1c")
    for t in th_clow + th_chigh:
        ax.axhline(t, ls="--", lw=0.8, color="#984ea3")
    ax.set_xlabel("Genes per cell")
    ax.set_ylabel("Total counts")
    ax.set_yscale("log")
    ax.set_title(f"Joint QC — {sample_id}")
    fig.tight_layout()
    fig.savefig(str(out.joint_png), dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    # ── 5. Summary TSV: ranges + per-threshold counts + percentile context ─
    feat_defs = [("n_genes_by_counts", ng, th_genes, "below"),
                 ("total_counts", tc, th_clow, "below"),
                 ("total_counts", tc, th_chigh, "above")]
    if mito_available:
        feat_defs.append(("pct_counts_mt", mt, th_mt, "above"))
    srows = []
    for col, vals, ths, op in feat_defs:
        fmin, fp50, fmax = float(np.min(vals)), float(np.median(vals)), float(np.max(vals))
        for t in ths:
            mask = (vals < t) if op == "below" else (vals > t)
            srows.append({"feature": col, "feature_min": round(fmin, 3),
                          "feature_p50": round(fp50, 3), "feature_max": round(fmax, 3),
                          "direction": op, "threshold": t,
                          "threshold_percentile": round(float((vals < t).mean() * 100), 3),
                          "n_removed": int(mask.sum()),
                          "pct_removed": round(100 * float(mask.mean()), 3),
                          "n_retained": int((~mask).sum())})
    pd.DataFrame(srows).to_csv(str(out.summary_tsv), sep="\t", index=False)

    # ── 6. Ingest-removed TSV (cell type × threshold) ────────────────────
    if ingest_vals is not None:
        irows = []
        for col, vals, ths, op in feat_defs:
            for t in ths:
                mask = (vals < t) if op == "below" else (vals > t)
                for ct, cnt in ingest_vals[mask].value_counts().items():
                    irows.append({"feature": col, "direction": op, "threshold": t,
                                  "ingest_celltype": ct, "n_removed": int(cnt)})
        pd.DataFrame(irows, columns=["feature", "direction", "threshold",
                                     "ingest_celltype", "n_removed"]
                     ).to_csv(str(out.ingest_tsv), sep="\t", index=False)
    else:
        with open(str(out.ingest_tsv), "w") as fh:
            fh.write("# no ingest reference supplied (qc_sweep.ingest_enabled=false)\n")
            fh.write("feature\tdirection\tthreshold\tingest_celltype\tn_removed\n")

    log.info("QC sweep complete for %s (no h5ad written by design)", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s",
              getattr(snakemake.params, "sample_id", "?"), traceback.format_exc())
    raise
