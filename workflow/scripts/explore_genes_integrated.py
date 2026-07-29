"""
explore_genes_integrated.py – Gene/signature exploration (integrated)
=====================================================================
Phase 1: loads the Harmony-integrated h5ad once, computes AUCell scores
for all signatures, writes expression_ranges.tsv, and produces PNGs
organised as:
    gene_exploration/{entry}/{entry}_genes.tsv
    gene_exploration/{entry}/Integrated/*.png
"""

import gc
import logging
import os
import sys
import tempfile
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image

from composition_barplots import find_niche_column, build_niche_palette
from explore_genes_common import (
    read_tsv_to_dict, classify_entries, apply_annotation_palette, apply_region_palette,
    make_score_adata, annotate_ct_region_dotplot, create_annotation_legend,
    composite_vertical, _compact_legend,
)

# Large cell type × region dotplots can exceed PIL's default pixel limit
Image.MAX_IMAGE_PIXELS = None

# Rasterise scatter data layers (text/axes stay vector)
sc.settings.set_figure_params(vector_friendly=True)


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
log = logging.getLogger("explore_genes_integrated")


# ── Helpers ──────────────────────────────────────────────────────────────────

def generate_sample_contribution(plot_adata, var_name, title_name, annot_key,
                                 has_regions, sample_col, out_path, dpi):
    """Across ALL samples, screen each sample's contribution to one entry's
    expression/score. Composite of dotplots (single entry) grouped by:
    sample, sample × region, sample × cell type. Stacked via PIL."""
    if sample_col is None or sample_col not in plot_adata.obs.columns:
        return
    ad = plot_adata.copy()
    samp = ad.obs[sample_col].astype(str)
    ad.obs["_sample"] = samp
    groupbys = [("_sample", "by sample")]
    if has_regions and "region_annotation" in ad.obs.columns:
        ad.obs["_sample_region"] = samp + " | " + ad.obs["region_annotation"].astype(str)
        groupbys.append(("_sample_region", "by sample \u00d7 region"))
    if annot_key in ad.obs.columns:
        ad.obs["_sample_ct"] = samp + " | " + ad.obs[annot_key].astype(str)
        groupbys.append(("_sample_ct", "by sample \u00d7 cell type"))

    tmp_files = []
    for gb, sub in groupbys:
        try:
            n_groups = ad.obs[gb].nunique()
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
            sc.pl.dotplot(ad, var_names=[var_name], groupby=gb,
                          standard_scale="var",
                          figsize=(max(5, 0.25 * n_groups + 3),
                                   max(3, n_groups * 0.3)),
                          title=f"{title_name} \u2013 {sub}", show=False)
            plt.savefig(tmp, dpi=dpi, bbox_inches="tight")
            plt.close("all")
            tmp_files.append(tmp)
        except Exception as e:
            log.warning("    Contribution dotplot (%s) failed: %s", sub, e)
            plt.close("all")
    if tmp_files:
        composite_vertical(tmp_files, out_path, dpi=dpi)
    for f in tmp_files:
        try:
            os.unlink(f)
        except OSError:
            pass


def generate_dotplots(adata, var_names, annot_key, has_regions, out_dir,
                      prefix, dpi, annotation_colors=None, region_colors=None):
    """Generate a composite PNG of dotplots (cell type, region, cell type ×
    region + annotation bars + legend) stacked vertically and right-aligned."""
    n_genes = len(var_names)
    fig_w = max(8, n_genes * 1.2 + 4)
    tmp_files = []

    # (1) By cell type
    try:
        n_ct = adata.obs[annot_key].nunique()
        tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
        sc.pl.dotplot(adata, var_names=var_names, groupby=annot_key,
                      standard_scale="var",
                      figsize=(fig_w, max(4, n_ct * 0.4)),
                      title=f"{prefix} – by cell type",
                      show=False)
        plt.savefig(tmp, dpi=dpi, bbox_inches="tight")
        plt.close("all")
        tmp_files.append(tmp)
    except Exception as e:
        log.warning("  Dotplot celltype failed: %s", e)
        plt.close("all")

    # (2) By region
    if has_regions:
        try:
            n_reg = adata.obs["region_annotation"].nunique()
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
            sc.pl.dotplot(adata, var_names=var_names,
                          groupby="region_annotation",
                          standard_scale="var",
                          figsize=(fig_w, max(4, n_reg * 0.5)),
                          title=f"{prefix} – by region",
                          show=False)
            plt.savefig(tmp, dpi=dpi, bbox_inches="tight")
            plt.close("all")
            tmp_files.append(tmp)
        except Exception as e:
            log.warning("  Dotplot region failed: %s", e)
            plt.close("all")

    # (3) By cell type × region (with row annotations)
    if has_regions:
        try:
            combined = (adata.obs[annot_key].astype(str) + " | "
                        + adata.obs["region_annotation"].astype(str))
            adata.obs["_ct_region"] = pd.Categorical(combined)
            n_groups = adata.obs["_ct_region"].nunique()
            fig_h = max(6, n_groups * 0.35)
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
            sc.pl.dotplot(adata, var_names=var_names, groupby="_ct_region",
                          standard_scale="var",
                          figsize=(fig_w, fig_h),
                          title=f"{prefix} – by cell type × region",
                          show=False)
            if annotation_colors or region_colors:
                annotate_ct_region_dotplot(annot_key,
                                          annotation_colors or {},
                                          region_colors or {})
            plt.savefig(tmp, dpi=dpi, bbox_inches="tight")
            tmp_files.append(tmp)

            # Standalone annotation legend (before plt.close so gcf() works)
            if annotation_colors or region_colors:
                tmp_leg = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
                create_annotation_legend(annot_key,
                                         annotation_colors or {},
                                         region_colors or {},
                                         tmp_leg, dpi)
                if os.path.isfile(tmp_leg) and os.path.getsize(tmp_leg) > 0:
                    tmp_files.append(tmp_leg)
            plt.close("all")
        except Exception as e:
            log.warning("  Dotplot celltype×region failed: %s", e)
            plt.close("all")
        finally:
            if "_ct_region" in adata.obs.columns:
                del adata.obs["_ct_region"]

    # Composite all dotplots + legend into a single PNG
    if tmp_files:
        composite_vertical(tmp_files,
                           os.path.join(out_dir, f"{prefix}_dotplots.png"),
                           dpi=dpi)
    for f in tmp_files:
        try:
            os.unlink(f)
        except OSError:
            pass


def generate_umap_composite(adata, color_col, annot_key, has_regions, out_path,
                             vmin, vmax, dpi, title="expression",
                             niche_col=None, niche_palette=None):
    """Composite PNG with 2–4 UMAPs (expression/score + cell type + region +
    spatial niche). The niche panel uses the shared palette and a compacted
    multi-column legend for high niche counts."""
    has_niche = bool(niche_col and niche_col in adata.obs.columns)
    n_panels = 2 + (1 if has_regions else 0) + (1 if has_niche else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 7),
                              gridspec_kw={"wspace": 0.5})

    try:
        sc.pl.umap(adata, color=color_col, size=2, frameon=False,
                    vmin=vmin, vmax=vmax, cmap="viridis",
                    title=title, ax=axes[0], show=False)
    except Exception as e:
        log.warning("  UMAP %s failed: %s", title, e)

    try:
        sc.pl.umap(adata, color=annot_key, size=2, frameon=False,
                    title="Cell types", legend_fontsize=6,
                    na_in_legend=False,
                    ax=axes[1], show=False)
    except Exception as e:
        log.warning("  UMAP celltype failed: %s", e)

    idx = 2
    if has_regions:
        try:
            sc.pl.umap(adata, color="region_annotation", size=2, frameon=False,
                        title="Regions", legend_fontsize=6,
                        na_in_legend=False,
                        ax=axes[idx], show=False)
        except Exception as e:
            log.warning("  UMAP region failed: %s", e)
        idx += 1

    if has_niche:
        try:
            if not isinstance(adata.obs[niche_col].dtype, pd.CategoricalDtype):
                adata.obs[niche_col] = adata.obs[niche_col].astype(str).astype("category")
            cats = list(adata.obs[niche_col].cat.categories)
            if niche_palette:
                adata.uns[f"{niche_col}_colors"] = [
                    niche_palette.get(str(c), "#cccccc") for c in cats]
            sc.pl.umap(adata, color=niche_col, size=2, frameon=False,
                        title="Spatial niches", legend_fontsize=5,
                        na_in_legend=False, ax=axes[idx], show=False)
            _compact_legend(axes[idx], title="Niche")
        except Exception as e:
            log.warning("  UMAP niche failed: %s", e)
        idx += 1

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ── Parameters ───────────────────────────────────────────────────────────────
integrated_path   = str(snakemake.input.integrated)
queries_path      = str(snakemake.input.queries)
out_ranges        = str(snakemake.output.ranges)
base_dir          = str(snakemake.params.outdir)

ANNOT_KEY         = str(snakemake.params.annot_key)
AUCELL_FRACTION   = float(snakemake.params.aucell_fraction)
NICHE_COLUMN      = str(snakemake.params.niche_column) if snakemake.params.niche_column else ""
DPI               = int(snakemake.params.dpi)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors
EXTRA_ANNOT_COLUMNS    = list(getattr(snakemake.params, "extra_annot_columns", []) or [])
SAMPLE_COLORS     = getattr(snakemake.params, "sample_colors", {}) or {}

try:
    log.info("=" * 70)
    log.info("Gene exploration – integrated plots")
    log.info("=" * 70)

    # ── 1. Parse queries ─────────────────────────────────────────────────
    entries_dict = read_tsv_to_dict(queries_path)
    individual_genes, signatures = classify_entries(entries_dict)
    all_entries = {**individual_genes, **signatures}
    log.info("Parsed %d entries: %d genes, %d signatures",
             len(all_entries), len(individual_genes), len(signatures))

    if not all_entries:
        log.warning("No entries — writing empty ranges and exiting.")
        pd.DataFrame(columns=["entry", "vmin", "vmax", "type"]).to_csv(
            out_ranges, sep="\t", index=False)
        sys.exit(0)

    # ── 2. Load integrated h5ad ──────────────────────────────────────────
    log.info("Loading integrated h5ad …")
    adata = sc.read_h5ad(integrated_path)
    log.info("  %d cells, %d genes", adata.n_obs, adata.n_vars)

    # Integrated UMAPs by each design column (grey for palette-less values)
    for _dc in EXTRA_ANNOT_COLUMNS:
        if _dc in adata.obs.columns:
            try:
                adata.obs[_dc] = adata.obs[_dc].astype(str).astype("category")
                cats = adata.obs[_dc].cat.categories
                pal = SAMPLE_COLORS.get(_dc, {}) if isinstance(SAMPLE_COLORS, dict) else {}
                adata.uns[f"{_dc}_colors"] = [pal.get(str(c), "#cccccc") for c in cats]
                _dd = os.path.join(base_dir, "_design")
                os.makedirs(_dd, exist_ok=True)
                sc.pl.umap(adata, color=[_dc], size=2, frameon=False,
                           title=f"Integrated – by {_dc}", show=False)
                plt.savefig(os.path.join(_dd, f"UMAP_by_{_dc}.png"),
                            dpi=DPI, bbox_inches="tight")
                plt.close()
            except Exception as e:
                log.warning("  design UMAP %s failed: %s", _dc, e)

    # ── 3. Validate genes ────────────────────────────────────────────────
    available_genes = set(adata.var_names)
    valid_entries = {}
    for name, gene_list in all_entries.items():
        present = [g for g in gene_list if g in available_genes]
        missing = [g for g in gene_list if g not in available_genes]
        if missing:
            log.warning("  %s: %d/%d genes missing: %s",
                        name, len(missing), len(gene_list), ", ".join(missing))
        if present:
            valid_entries[name] = present
        else:
            log.warning("  %s: ALL genes missing — skipping.", name)
    individual_genes = {k: v for k, v in individual_genes.items() if k in valid_entries}
    signatures = {k: v for k, v in signatures.items() if k in valid_entries}
    log.info("  %d valid entries", len(valid_entries))

    # ── 4. AUCell for signatures ─────────────────────────────────────────
    aucell_scores = {}
    if signatures:
        log.info("Computing AUCell for %d signatures …", len(signatures))
        import decoupler as dc

        net_rows = []
        for sig_name in signatures:
            for gene in valid_entries[sig_name]:
                net_rows.append({"source": sig_name, "target": gene, "weight": 1.0})
        net_df = pd.DataFrame(net_rows)
        n_up = max(1, int(adata.n_vars * AUCELL_FRACTION))
        log.info("  n_up = %d (%.1f%% of %d genes)",
                 n_up, AUCELL_FRACTION * 100, adata.n_vars)

        dc.run_aucell(adata, net=net_df, source="source", target="target",
                      n_up=n_up, use_raw=False, verbose=True)

        if "aucell_estimate" in adata.obsm:
            est = adata.obsm["aucell_estimate"]
            for sig_name in signatures:
                if sig_name in est.columns:
                    aucell_scores[sig_name] = est[sig_name].values
                    adata.obs[f"AUCell_{sig_name}"] = est[sig_name].values

    # ── 5. Expression ranges ─────────────────────────────────────────────
    log.info("Computing expression ranges …")
    ranges_rows = []

    # Collect ALL unique genes (standalone + signature members)
    all_unique_genes = set(individual_genes.keys())
    for sig_name in signatures:
        for gene in valid_entries[sig_name]:
            all_unique_genes.add(gene)

    for gene_name in sorted(all_unique_genes):
        if gene_name not in adata.var_names:
            log.warning("  gene '%s' not in integrated object — skipping range", gene_name)
            continue
        expr = adata[:, gene_name].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        expr = np.asarray(expr).flatten()
        vmin, vmax = float(np.percentile(expr, 1)), float(np.percentile(expr, 99))
        ranges_rows.append({"entry": gene_name, "vmin": vmin, "vmax": vmax, "type": "gene"})

    for sig_name in signatures:
        if sig_name in aucell_scores:
            scores = aucell_scores[sig_name]
            vmin = float(np.percentile(scores, 1))
            vmax = float(np.percentile(scores, 99))
            ranges_rows.append({"entry": sig_name, "vmin": vmin, "vmax": vmax,
                                "type": "signature"})

    ranges_df = pd.DataFrame(ranges_rows)
    ranges_df.to_csv(out_ranges, sep="\t", index=False)
    log.info("  Written → %s", out_ranges)

    # ── 6. Apply palettes ────────────────────────────────────────────────
    apply_annotation_palette(adata, ANNOT_KEY, ANNOTATION_COLORS)
    apply_region_palette(adata, REGION_COLORS)
    if ANNOT_KEY in adata.obs.columns:
        adata.obs[ANNOT_KEY] = adata.obs[ANNOT_KEY].astype("category")

    has_regions = ("region_annotation" in adata.obs.columns
                   and adata.obs["region_annotation"].nunique() > 1
                   and not all(adata.obs["region_annotation"] == "Unlabeled"))
    has_niche = bool(NICHE_COLUMN) and NICHE_COLUMN in adata.obs.columns
    sample_col = next((c for c in ["sample", "sample_batch"]
                       if c in adata.obs.columns), None)

    # Spatial-niche palette (shared, value-deterministic) → niche colours match
    # the spatial_niches rule in the UMAP composite niche panel.
    niche_col = NICHE_COLUMN if has_niche else None
    niche_palette = None
    if has_niche:
        _nv = adata.obs[niche_col].astype(str)
        try:
            _ncats = sorted(_nv.unique(), key=lambda x: int(x))
        except ValueError:
            _ncats = sorted(_nv.unique())
        niche_palette = build_niche_palette(
            _ncats, ANNOTATION_COLORS.get("spatial_niche", {})
            if isinstance(ANNOTATION_COLORS, dict) else {})
        log.info("spatial niche '%s' (%d niches) → UMAP panel + dotplots",
                 niche_col, len(_ncats))

    # ── 7. Generate PNGs ─────────────────────────────────────────────────

    for gene_name in individual_genes:
        log.info("Plotting gene: %s", gene_name)
        entry_dir = os.path.join(base_dir, gene_name)
        int_dir = os.path.join(entry_dir, "Integrated")
        os.makedirs(int_dir, exist_ok=True)

        # Provenance TSV
        pd.DataFrame({gene_name: [gene_name]}).to_csv(
            os.path.join(entry_dir, f"{gene_name}_genes.tsv"),
            sep="\t", index=False)

        row = ranges_df[ranges_df["entry"] == gene_name].iloc[0]
        vmin, vmax = row["vmin"], row["vmax"]

        generate_dotplots(adata, [gene_name], ANNOT_KEY, has_regions,
                          int_dir, gene_name, DPI, ANNOTATION_COLORS, REGION_COLORS)
        generate_umap_composite(adata, gene_name, ANNOT_KEY, has_regions,
                                os.path.join(int_dir, f"{gene_name}_UMAPs.png"),
                                vmin, vmax, DPI, "expression",
                                niche_col=niche_col, niche_palette=niche_palette)

        # Per-entry, all-samples contribution (sample / sample×area / sample×ct)
        generate_sample_contribution(
            adata, gene_name, gene_name, ANNOT_KEY, has_regions, sample_col,
            os.path.join(int_dir, f"{gene_name}_sample_contribution.png"), DPI)

        if has_niche:
            try:
                n_niche = adata.obs[NICHE_COLUMN].nunique()
                sc.pl.dotplot(adata, var_names=[gene_name],
                              groupby=NICHE_COLUMN, standard_scale="var",
                              figsize=(8, max(4, n_niche * 0.5)),
                              title=f"{gene_name} – by niche", show=False)
                plt.savefig(os.path.join(int_dir,
                            f"{gene_name}_dotplot_niche.png"),
                            dpi=DPI, bbox_inches="tight")
                plt.close("all")
            except Exception as e:
                log.warning("  Niche dotplot failed: %s", e)
                plt.close("all")

    for sig_name in signatures:
        log.info("Plotting signature: %s", sig_name)
        score_col = f"AUCell_{sig_name}"
        if score_col not in adata.obs.columns:
            log.warning("  AUCell score not found — skipping.")
            continue

        entry_dir = os.path.join(base_dir, sig_name)
        int_dir = os.path.join(entry_dir, "Integrated")
        os.makedirs(int_dir, exist_ok=True)

        present_genes = valid_entries[sig_name]
        pd.DataFrame({sig_name: present_genes}).to_csv(
            os.path.join(entry_dir, f"{sig_name}_genes.tsv"),
            sep="\t", index=False)

        row = ranges_df[ranges_df["entry"] == sig_name].iloc[0]
        vmin, vmax = row["vmin"], row["vmax"]

        # Member-gene dotplots
        generate_dotplots(adata, present_genes, ANNOT_KEY, has_regions,
                          int_dir, f"{sig_name}_genes", DPI,
                          ANNOTATION_COLORS, REGION_COLORS)

        # AUCell score dotplots
        score_ad = make_score_adata(adata, score_col)
        generate_dotplots(score_ad, [score_col], ANNOT_KEY, has_regions,
                          int_dir, f"{sig_name}_aucell", DPI,
                          ANNOTATION_COLORS, REGION_COLORS)

        # Per-entry, all-samples contribution of the AUCell score
        generate_sample_contribution(
            score_ad, score_col, f"{sig_name} AUCell", ANNOT_KEY, has_regions,
            sample_col,
            os.path.join(int_dir, f"{sig_name}_sample_contribution.png"), DPI)
        del score_ad

        # UMAPs
        generate_umap_composite(adata, score_col, ANNOT_KEY, has_regions,
                                os.path.join(int_dir, f"{sig_name}_UMAPs.png"),
                                vmin, vmax, DPI, "AUCell",
                                niche_col=niche_col, niche_palette=niche_palette)

        if has_niche:
            try:
                score_ad = make_score_adata(adata, score_col)
                n_niche = adata.obs[NICHE_COLUMN].nunique()
                sc.pl.dotplot(score_ad, var_names=[score_col],
                              groupby=NICHE_COLUMN, standard_scale="var",
                              figsize=(8, max(4, n_niche * 0.5)),
                              title=f"{sig_name} AUCell – by niche",
                              show=False)
                plt.savefig(os.path.join(int_dir,
                            f"{sig_name}_aucell_dotplot_niche.png"),
                            dpi=DPI, bbox_inches="tight")
                plt.close("all")
                del score_ad
            except Exception as e:
                log.warning("  Niche dotplot failed: %s", e)
                plt.close("all")

        # Member-gene expression per niche (genes × niches) — new separate plot
        if has_niche and present_genes:
            try:
                pg = [g for g in present_genes if g in adata.var_names]
                if pg:
                    n_niche = adata.obs[NICHE_COLUMN].nunique()
                    fig_w = max(8, len(pg) * 0.5 + 4)
                    sc.pl.dotplot(adata, var_names=pg, groupby=NICHE_COLUMN,
                                  standard_scale="var",
                                  figsize=(fig_w, max(4, n_niche * 0.5)),
                                  title=f"{sig_name} genes – by niche",
                                  show=False)
                    plt.savefig(os.path.join(int_dir,
                                f"{sig_name}_genes_dotplot_niche.png"),
                                dpi=DPI, bbox_inches="tight")
                    plt.close("all")
            except Exception as e:
                log.warning("  Niche member-gene dotplot failed: %s", e)
                plt.close("all")

    # ── Overall AUCell-by-niche overview: ALL signatures' scores × niches in
    #    one dotplot (overall contribution of each signature per niche) ──────
    if has_niche and signatures:
        try:
            sig_cols = [f"AUCell_{s}" for s in signatures
                        if f"AUCell_{s}" in adata.obs.columns]
            if sig_cols:
                score_ad = sc.AnnData(
                    X=adata.obs[sig_cols].to_numpy(dtype=np.float32),
                    obs=adata.obs.drop(columns=sig_cols, errors="ignore"),
                    var=pd.DataFrame(
                        index=[c.replace("AUCell_", "") for c in sig_cols]))
                n_niche = adata.obs[NICHE_COLUMN].nunique()
                fig_w = max(8, len(sig_cols) * 0.5 + 4)
                sc.pl.dotplot(score_ad, var_names=list(score_ad.var_names),
                              groupby=NICHE_COLUMN, standard_scale="var",
                              figsize=(fig_w, max(4, n_niche * 0.5)),
                              title="AUCell signatures – by niche (overall)",
                              show=False)
                plt.savefig(os.path.join(base_dir, "aucell_overall_by_niche.png"),
                            dpi=DPI, bbox_inches="tight")
                plt.close("all")
                del score_ad
        except Exception as e:
            log.warning("Overall AUCell-by-niche dotplot failed: %s", e)
            plt.close("all")

    del adata
    gc.collect()
    log.info("Done.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
