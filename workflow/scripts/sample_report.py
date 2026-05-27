"""
sample_report.py – Multi-page PDF report, one page per sample
================================================================
Each page contains:
  Row 1:  UMAP (clusters) | UMAP (annotation) | Spatial (clusters) | Spatial (annotation)
  Row 2:  Dotplot (cell type markers)          | Barplot (cell proportions by region)

Memory is bounded by loading one sample at a time.
"""

import csv
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

def read_tsv_to_dict(tsv_path):
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = {field: [] for field in reader.fieldnames}
        for row in reader:
            for key, val in row.items():
                if val.strip():
                    columns[key].append(val.strip())
    return columns


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


# ── Parameters ───────────────────────────────────────────────────────────────
annotated_paths   = [str(p) for p in snakemake.input.annotated]
markers_path      = str(snakemake.input.annotation_markers)
sample_ids        = list(snakemake.params.sample_ids)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors
out_report        = str(snakemake.output.report)

try:
    log.info("=" * 70)
    log.info("Generating sample report PDF for %d samples", len(sample_ids))
    log.info("=" * 70)

    Path(out_report).parent.mkdir(parents=True, exist_ok=True)

    # Read marker genes
    marker_dict = {}
    if markers_path and os.path.isfile(markers_path):
        raw_markers = read_tsv_to_dict(markers_path)
        marker_dict = {ct: genes for ct, genes in raw_markers.items() if genes}

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

                # Find best leiden column
                leiden_cols = sorted([c for c in adata.obs.columns if c.startswith("leiden_")])
                leiden_key = leiden_cols[0] if leiden_cols else None

                # Find best annotation column
                annot_key = None
                for candidate in ["cell_type_tsv", "cell_type_refined",
                                  "cell_type_ingest", "cell_type_external"]:
                    if candidate in adata.obs.columns:
                        annot_key = candidate
                        break

                # Apply palettes
                for obs_key in ["cell_type_tsv", "cell_type_refined",
                                "cell_type_ingest", "cell_type_external"]:
                    apply_palette(adata, obs_key, ANNOTATION_COLORS)
                apply_region_palette(adata, REGION_COLORS)

                # ── Build figure ─────────────────────────────────────────
                fig = plt.figure(figsize=(28, 16))
                fig.suptitle(f"Sample: {sample_id}", fontsize=18, fontweight="bold", y=0.98)
                gs = gridspec.GridSpec(2, 4, figure=fig,
                                       height_ratios=[1, 1.2],
                                       hspace=0.35, wspace=0.35)

                # Row 1: UMAP (clusters) | UMAP (annotation) | Spatial (clusters) | Spatial (annotation)
                # Panel 1: UMAP clusters
                if leiden_key:
                    ax1 = fig.add_subplot(gs[0, 0])
                    try:
                        sc.pl.umap(adata, color=leiden_key, size=2, frameon=False,
                                   title="Clusters", ax=ax1, show=False)
                    except Exception as e:
                        log.warning("  UMAP clusters failed: %s", e)
                        ax1.set_title("Clusters (failed)")

                # Panel 2: UMAP annotation
                if annot_key:
                    ax2 = fig.add_subplot(gs[0, 1])
                    try:
                        sc.pl.umap(adata, color=annot_key, size=2, frameon=False,
                                   title=annot_key.replace("cell_type_", ""),
                                   ax=ax2, show=False)
                    except Exception as e:
                        log.warning("  UMAP annotation failed: %s", e)
                        ax2.set_title("Annotation (failed)")

                # Panel 3: Spatial clusters
                if leiden_key:
                    ax3 = fig.add_subplot(gs[0, 2])
                    try:
                        sc.pl.spatial(adata, color=leiden_key, spot_size=20, frameon=False,
                                      title="Spatial – clusters",
                                      library_id=library_id, ax=ax3, show=False)
                    except Exception as e:
                        log.warning("  Spatial clusters failed: %s", e)
                        ax3.set_title("Spatial clusters (failed)")

                # Panel 4: Spatial annotation
                if annot_key:
                    ax4 = fig.add_subplot(gs[0, 3])
                    try:
                        sc.pl.spatial(adata, color=annot_key, spot_size=20, frameon=False,
                                      title="Spatial – " + annot_key.replace("cell_type_", ""),
                                      library_id=library_id, ax=ax4, show=False)
                    except Exception as e:
                        log.warning("  Spatial annotation failed: %s", e)
                        ax4.set_title("Spatial annotation (failed)")

                # Row 2: Dotplot (spanning 3 columns) | Barplot (1 column)

                # Panel 5: Dotplot with cell type markers
                ax5 = fig.add_subplot(gs[1, 0:3])
                if annot_key and marker_dict:
                    try:
                        # Filter markers to genes present in adata
                        filtered = {}
                        for ct, genes in marker_dict.items():
                            present = [g for g in genes if g in adata.var_names]
                            if present:
                                filtered[ct] = present

                        if filtered:
                            # Ensure annotation is category with unused removed
                            adata.obs[annot_key] = adata.obs[annot_key].cat.remove_unused_categories()
                            sc.pl.dotplot(
                                adata, filtered, groupby=annot_key,
                                standard_scale="var", swap_axes=True,
                                title="Cell type markers", ax=ax5, show=False)
                        else:
                            ax5.text(0.5, 0.5, "No marker genes found in adata",
                                     ha="center", va="center", transform=ax5.transAxes)
                            ax5.set_title("Cell type markers")
                    except Exception as e:
                        log.warning("  Dotplot failed: %s", e)
                        ax5.text(0.5, 0.5, f"Dotplot failed: {e}",
                                 ha="center", va="center", transform=ax5.transAxes,
                                 fontsize=8)
                else:
                    ax5.text(0.5, 0.5, "No annotation or markers available",
                             ha="center", va="center", transform=ax5.transAxes)
                    ax5.set_title("Cell type markers")

                # Panel 6: Composition barplot (by region)
                ax6 = fig.add_subplot(gs[1, 3])
                has_regions = ("region_annotation" in adata.obs.columns
                               and adata.obs["region_annotation"].nunique() > 1
                               and not all(adata.obs["region_annotation"] == "Unlabeled"))
                if annot_key and has_regions:
                    try:
                        ct = pd.crosstab(adata.obs["region_annotation"],
                                         adata.obs[annot_key])
                        ct_norm = ct.div(ct.sum(axis=1), axis=0) * 100

                        colors = None
                        if isinstance(ANNOTATION_COLORS, dict):
                            cd = ANNOTATION_COLORS.get(annot_key, {})
                            if cd:
                                colors = [cd.get(str(c), "#cccccc") for c in ct_norm.columns]

                        ct_norm.plot(kind="bar", stacked=True, ax=ax6, color=colors,
                                     legend=False)
                        ax6.set_ylabel("Percentage (%)")
                        ax6.set_xlabel("Region")
                        ax6.set_title("Cell proportions by region")
                        ax6.tick_params(axis="x", rotation=45)
                        # Compact legend below the barplot
                        ax6.legend(bbox_to_anchor=(1.02, 1), loc="upper left",
                                   fontsize=6, title=annot_key.replace("cell_type_", ""),
                                   title_fontsize=7)
                    except Exception as e:
                        log.warning("  Barplot failed: %s", e)
                        ax6.text(0.5, 0.5, f"Barplot failed: {e}",
                                 ha="center", va="center", transform=ax6.transAxes,
                                 fontsize=8)
                else:
                    ax6.text(0.5, 0.5, "No region annotation",
                             ha="center", va="center", transform=ax6.transAxes)
                    ax6.set_title("Cell proportions by region")

                pdf.savefig(fig, dpi=200, bbox_inches="tight")
                plt.close(fig)
                log.info("  Page saved for %s", sample_id)

            except Exception as e:
                log.warning("FAILED page for %s: %s\n%s", sample_id, e,
                            traceback.format_exc())
                # Create a placeholder page with error message
                fig, ax = plt.subplots(figsize=(28, 16))
                ax.text(0.5, 0.5, f"Sample {sample_id}\nFailed: {e}",
                        ha="center", va="center", fontsize=14)
                ax.set_axis_off()
                pdf.savefig(fig, dpi=200)
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
