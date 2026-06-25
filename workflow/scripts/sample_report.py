"""
sample_report.py – Multi-page PDF report, one page per sample
================================================================
Each page contains:
  Row 1:  UMAP (clusters) | UMAP (annotation) | Spatial (clusters) | Spatial (annotation)
  Row 2:  Dotplot (cell type markers)          | Barplot (cell proportions by region)

Memory is bounded by loading one sample at a time.
"""

import gc
import logging
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scanpy as sc

# Shared composition-barplot helpers (black edge, legend==stack order, config colours)
try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:                      # very old Snakemake
    _here = os.getcwd()
sys.path.insert(0, _here)
from composition_barplots import (draw_stacked_composition, find_niche_column,
                                  build_niche_palette)


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
log = logging.getLogger("sample_report")


# ── Helpers ──────────────────────────────────────────────────────────────────

def _rasterize_heavy_layers(fig):
    """Rasterise the heavy scatter layers (UMAP / spatial spot collections and
    any embedded images) of every axis in a figure, keeping text, lines, axes
    and legends as vector. This works around scanpy's ``vector_friendly``
    rasterisation being applied inconsistently (scverse/scanpy#2005), where
    spatial scatters in particular stay vector and bloat the PDF (one path per
    spot → hundreds of MB). When the page is then written with a finite dpi, the
    rasterised layers are baked at that resolution, cutting file size by ~10x
    while labels stay crisp/selectable."""
    import matplotlib.collections as mcoll
    import matplotlib.image as mimage
    for ax in fig.get_axes():
        for art in ax.get_children():
            if isinstance(art, (mcoll.Collection, mimage.AxesImage)):
                art.set_rasterized(True)


def apply_palette(adata, obs_key, colors_cfg):
    """Apply annotation_colors palette to an obs column."""
    if not isinstance(colors_cfg, dict):
        return
    cd = colors_cfg.get(obs_key, {})
    if cd and obs_key in adata.obs.columns:
        adata.obs[obs_key] = adata.obs[obs_key].astype("category")
        cats = adata.obs[obs_key].cat.categories
        adata.uns[f"{obs_key}_colors"] = [cd.get(str(c), "#cccccc") for c in cats]


def apply_region_palette(adata, region_colors):
    """Apply region_colors palette to region_annotation."""
    if region_colors and "region_annotation" in adata.obs.columns:
        adata.obs["region_annotation"] = adata.obs["region_annotation"].astype("category")
        cats = adata.obs["region_annotation"].cat.categories
        adata.uns["region_annotation_colors"] = [
            region_colors.get(str(c), "#cccccc") for c in cats
        ]


def _pick_keys(adata, niche_column=None):
    """Return (manual_leiden, tsv, auto, niche) obs keys (any may be None)."""
    manual = adata.uns.get("annotation_leiden_key")
    if manual not in adata.obs.columns:
        leiden_cols = sorted(c for c in adata.obs.columns if c.startswith("leiden_"))
        manual = leiden_cols[0] if leiden_cols else None
    tsv = "cell_type_tsv" if "cell_type_tsv" in adata.obs.columns else None
    auto = next((c for c in ["cell_type_ingest", "cell_type_external"]
                 if c in adata.obs.columns), None)
    niche = find_niche_column(adata, niche_column)
    return manual, tsv, auto, niche


def _build_page1(adata, sample_id, keys, library_id):
    """Row 1 = UMAPs, Row 2 = spatial, for clusters / tsv / auto / niche."""
    manual, tsv, auto, niche = keys
    panels = [(manual, "Clusters"), (tsv, "TSV annotation"),
              (auto, "Auto annotation"), (niche, "Spatial niche")]
    panels = [(k, t) for k, t in panels if k]
    n = len(panels)
    fig, axes = plt.subplots(2, n, figsize=(8 * n, 14),
                             gridspec_kw={"wspace": 0.5, "hspace": 0.25},
                             squeeze=False)
    fig.suptitle(f"Sample: {sample_id}", fontsize=18, fontweight="bold", y=0.99)
    for j, (key, title) in enumerate(panels):
        try:
            sc.pl.umap(adata, color=key, size=2, frameon=False, title=title,
                       legend_fontsize=6, na_in_legend=False,
                       ax=axes[0, j], show=False)
        except Exception as e:
            log.warning("  UMAP %s failed: %s", title, e); axes[0, j].set_title(f"{title} (failed)")
        try:
            sc.pl.spatial(adata, color=key, spot_size=20, frameon=False,
                          title=title, library_id=library_id,
                          legend_fontsize=6, na_in_legend=False,
                          ax=axes[1, j], show=False)
        except Exception as e:
            log.warning("  Spatial %s failed: %s", title, e); axes[1, j].set_title(f"{title} (failed)")
    return fig


def _top_markers_by_tsv(adata, tsv_key, n=10):
    """Top-n DE markers per tsv cell type (deduped, order preserved)."""
    try:
        adata.obs[tsv_key] = adata.obs[tsv_key].astype("category")
        adata.obs[tsv_key] = adata.obs[tsv_key].cat.remove_unused_categories()
        if adata.obs[tsv_key].nunique() < 2:
            return []
        sc.tl.rank_genes_groups(adata, groupby=tsv_key, method="wilcoxon")
        names = adata.uns["rank_genes_groups"]["names"]
        ordered = []
        for grp in names.dtype.names:
            for g in list(names[grp])[:n]:
                if g in adata.var_names and g not in ordered:
                    ordered.append(g)
        return ordered
    except Exception as e:
        log.warning("  marker ranking failed: %s", e)
        return []


def _build_page2(adata, sample_id, keys, annotation_colors, sample_col,
                 has_regions):
    """Return a list of page-2 figures.

    • Composition barplots: when the sample has region and/or spatial-niche
      information, the cell-type composition is broken down by those (one
      absolute + one relative panel each). Only when NEITHER is present does it
      fall back to a single by-sample breakdown (the otherwise-trivial view).
    • A horizontal marker dotplot (genes on X) with a compact vertical colour
      bar and an accurately-sized dot-size legend.
    """
    manual, tsv, auto, niche = keys
    cmap = (annotation_colors.get(tsv, {})
            if (tsv and isinstance(annotation_colors, dict)) else {})
    figs = []

    # Group the composition by region/niche when available; only fall back to
    # the (trivial) per-sample breakdown when there is no region and no niche.
    bases = []
    if has_regions:
        bases.append(("region_annotation", "Region"))
    if niche:
        bases.append((niche, "Niche"))
    if not bases:
        bases.append((sample_col, "Sample"))

    if tsv:
        try:
            ncol = 2 * len(bases)
            figb = plt.figure(figsize=(6.5 * ncol, 5.5))
            figb.suptitle(f"Sample: {sample_id} — composition",
                          fontsize=16, fontweight="bold", y=1.02)
            gs = gridspec.GridSpec(1, ncol, figure=figb, wspace=0.5)
            for i, (gkey, glabel) in enumerate(bases):
                ct = pd.crosstab(adata.obs[gkey], adata.obs[tsv])
                ax_abs = figb.add_subplot(gs[0, 2 * i])
                draw_stacked_composition(ax_abs, ct, cmap, normalize=False,
                                         ylabel="Number of cells", xlabel=glabel,
                                         title=f"by {glabel} (absolute)",
                                         legend=False)
                ax_rel = figb.add_subplot(gs[0, 2 * i + 1])
                draw_stacked_composition(ax_rel, ct, cmap, normalize=True,
                                         ylabel="Percentage (%)", xlabel=glabel,
                                         title=f"by {glabel} (relative)",
                                         legend_title="Cell type",
                                         legend=(i == len(bases) - 1))
            figs.append(figb)
        except Exception as e:
            log.warning("  composition barplots failed: %s", e)

    # Spatial-niche composition per region (niche AS the stack; colours match
    # the spatial_niches rule via the shared value-deterministic palette).
    if niche and has_regions:
        try:
            nvals = adata.obs[niche].astype(str)
            try:
                ncats = sorted(nvals.unique(), key=lambda x: int(x))
            except ValueError:
                ncats = sorted(nvals.unique())
            niche_cfg = (annotation_colors.get("spatial_niche", {})
                         if isinstance(annotation_colors, dict) else {})
            niche_pal = build_niche_palette(ncats, niche_cfg)
            ctn = pd.crosstab(adata.obs["region_annotation"], nvals)
            fign = plt.figure(figsize=(13, 5.5))
            fign.suptitle(f"Sample: {sample_id} — spatial niche composition",
                          fontsize=16, fontweight="bold", y=1.02)
            gsn = gridspec.GridSpec(1, 2, figure=fign, wspace=0.5)
            axn_abs = fign.add_subplot(gsn[0, 0])
            draw_stacked_composition(axn_abs, ctn, niche_pal, normalize=False,
                                     ylabel="Number of cells", xlabel="Region",
                                     title="Niche by region (absolute)",
                                     legend=False)
            axn_rel = fign.add_subplot(gsn[0, 1])
            draw_stacked_composition(axn_rel, ctn, niche_pal, normalize=True,
                                     ylabel="Percentage (%)", xlabel="Region",
                                     title="Niche by region (relative)",
                                     legend_title="Niche",
                                     legend_kwargs=dict(
                                         ncol=max(1, len(ncats) // 15 + 1),
                                         fontsize=6))
            figs.append(fign)
        except Exception as e:
            log.warning("  niche-by-region barplot failed: %s", e)

    # ── Horizontal marker dotplot (genes on X), native scanpy legend column ──
    genes = _top_markers_by_tsv(adata, tsv, n=10) if tsv else []
    if genes and tsv:
        try:
            n_g = len(genes)
            n_grp = adata.obs[tsv].astype("category").nunique()
            dp = sc.pl.dotplot(
                adata, var_names=genes, groupby=tsv, standard_scale="var",
                swap_axes=False,                          # genes on X → horizontal
                figsize=(max(9, n_g * 0.34), max(3.0, n_grp * 0.45)),
                title=f"Sample: {sample_id} — top-10 markers per TSV cell type",
                return_fig=True)
            dp.legend(width=1.2)                          # compact, vertical colourbar
            dp.make_figure()
            main_ax = dp.ax_dict.get("mainplot_ax")
            if main_ax is not None:
                for lbl in main_ax.get_xticklabels():
                    lbl.set_rotation(45)
                    lbl.set_ha("right")
            figs.append(dp.fig)
        except Exception as e:
            log.warning("  marker dotplot failed: %s", e)
            ferr, axe = plt.subplots(figsize=(11, 4))
            axe.text(0.5, 0.5, f"Dotplot failed: {e}", ha="center", va="center",
                     transform=axe.transAxes, fontsize=8)
            axe.set_axis_off()
            figs.append(ferr)

    if not figs:
        f, ax = plt.subplots(figsize=(11, 6))
        ax.text(0.5, 0.5, "No markers / composition to display",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        figs.append(f)
    return figs


# ── Parameters ───────────────────────────────────────────────────────────────
annotated_paths   = [str(p) for p in snakemake.input.annotated]
sample_ids        = list(snakemake.params.sample_ids)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors
DPI          = int(getattr(snakemake.params, "dpi", 300))
NICHE_COLUMN      = getattr(snakemake.params, "niche_column", "")
out_report        = str(snakemake.output.report)

try:
    log.info("=" * 70)
    log.info("Generating sample report PDF for %d samples", len(sample_ids))
    log.info("=" * 70)

    Path(out_report).parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(out_report) as pdf:
        for adata_path, sample_id in zip(annotated_paths, sample_ids):
            log.info("Processing sample: %s", sample_id)
            try:
                adata = sc.read_h5ad(adata_path)
                log.info("  Loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)

                # Determine library_id for spatial plots
                library_id = None
                if "spatial" in adata.uns and len(adata.uns["spatial"]) == 1:
                    library_id = list(adata.uns["spatial"].keys())[0]

                # palettes
                for obs_key in ["cell_type_tsv", "cell_type_ingest",
                                "cell_type_external"]:
                    apply_palette(adata, obs_key, ANNOTATION_COLORS)
                apply_region_palette(adata, REGION_COLORS)

                sample_col = next((c for c in ["sample", "sample_batch"]
                                   if c in adata.obs.columns), None)
                if sample_col is None:
                    adata.obs["sample"] = sample_id
                    sample_col = "sample"

                keys = _pick_keys(adata, NICHE_COLUMN)

                has_regions = ("region_annotation" in adata.obs.columns
                               and adata.obs["region_annotation"].nunique() > 1
                               and not all(adata.obs["region_annotation"] == "Unlabeled"))

                fig1 = _build_page1(adata, sample_id, keys, library_id)
                _rasterize_heavy_layers(fig1)
                pdf.savefig(fig1, bbox_inches="tight", dpi=DPI); plt.close(fig1)

                for fig2 in _build_page2(adata, sample_id, keys,
                                         ANNOTATION_COLORS, sample_col, has_regions):
                    _rasterize_heavy_layers(fig2)
                    pdf.savefig(fig2, bbox_inches="tight", dpi=DPI); plt.close(fig2)
                log.info("  Pages saved for %s", sample_id)

            except Exception as e:
                log.warning("FAILED page for %s: %s\n%s", sample_id, e,
                            traceback.format_exc())
                # Create a placeholder page with error message
                fig, ax = plt.subplots(figsize=(28, 16))
                ax.text(0.5, 0.5, f"Sample {sample_id}\nFailed: {e}",
                        ha="center", va="center", fontsize=14)
                ax.set_axis_off()
                pdf.savefig(fig, dpi=DPI)
                plt.close(fig)
            finally:
                # Free memory before next sample
                if "adata" in dir():
                    del adata
                gc.collect()

    log.info("Report saved → %s", out_report)

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
