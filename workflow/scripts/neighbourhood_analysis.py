"""
neighbourhood_analysis.py – Spatial neighbourhood analysis with squidpy
=========================================================================
Uses the cell_type column from the annotated adata based on annot_type.
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
import squidpy as sq


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
log = logging.getLogger("neighbourhood_analysis")

sample_id   = snakemake.params.sample_id
annot_type  = snakemake.params.annot_type
annotation_colors = dict(getattr(snakemake.params, "annotation_colors", {}) or {})
# Cell-type annotation columns inherit cell_type_tsv (coherence for external/ingest/refined).
if annotation_colors.get("cell_type_tsv"):
    for _k in ("cell_type_external", "cell_type_ingest", "cell_type_refined"):
        annotation_colors.setdefault(_k, annotation_colors["cell_type_tsv"])
adata_path  = str(snakemake.input.adata)
results_dir = str(snakemake.output.results_dir)

ANNOT_COL_MAP = {
    "tsv_annotation": "cell_type_tsv",
    "refined_annotation": "cell_type_refined",
    "external_annotation": "cell_type_external",
}

try:
    log.info("=" * 70)
    log.info("Neighbourhood: sample=%s, annot_type=%s", sample_id, annot_type)
    log.info("=" * 70)

    os.makedirs(results_dir, exist_ok=True)

    cell_type_col = ANNOT_COL_MAP.get(annot_type)
    if cell_type_col is None:
        raise ValueError(f"Unknown annot_type '{annot_type}'.")

    adata = sc.read_h5ad(adata_path)

    if cell_type_col not in adata.obs.columns:
        log.warning("Column '%s' not found. Skipping.", cell_type_col)
        Path(os.path.join(results_dir, f"SKIPPED_no_{cell_type_col}.txt")).write_text(
            f"Column {cell_type_col} not found.\n"
        )
        sys.exit(0)

    adata.obs["cell_type"] = adata.obs[cell_type_col].astype("category")
    adata = adata[adata.obs["cell_type"] != "Unannotated"].copy()

    if adata.n_obs < 50 or adata.obs["cell_type"].nunique() < 2:
        log.warning("Too few cells (%d) or types (%d). Skipping.",
                     adata.n_obs, adata.obs["cell_type"].nunique())
        Path(os.path.join(results_dir, "SKIPPED_insufficient_data.txt")).write_text(
            f"cells={adata.n_obs}, types={adata.obs['cell_type'].nunique()}\n"
        )
        sys.exit(0)

    log.info("%d cells, %d cell types", adata.n_obs, adata.obs["cell_type"].nunique())

    # Map the config palette onto cell_type so squidpy plots (co-occurrence lines,
    # enrichment axes) use consistent, publication colours. Grey fallback for
    # values absent from the palette.
    adata.obs["cell_type"] = adata.obs["cell_type"].cat.remove_unused_categories()
    _ct_pal_cfg = annotation_colors.get(cell_type_col, {}) if isinstance(annotation_colors, dict) else {}
    cell_type_palette = {str(c): _ct_pal_cfg.get(str(c), "#cccccc")
                         for c in adata.obs["cell_type"].cat.categories}
    adata.uns["cell_type_colors"] = [cell_type_palette[str(c)]
                                     for c in adata.obs["cell_type"].cat.categories]

    # Global neighbourhood analysis
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    sq.gr.nhood_enrichment(adata, cluster_key="cell_type")

    sq.pl.nhood_enrichment(adata, cluster_key="cell_type", figsize=(8, 8))
    plt.savefig(os.path.join(results_dir, f"nhood_enrichment_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    if "cell_type_nhood_enrichment" in adata.uns:
        zscore = adata.uns["cell_type_nhood_enrichment"]["zscore"]
        cats = adata.obs["cell_type"].cat.categories
        pd.DataFrame(zscore, index=cats, columns=cats).to_csv(
            os.path.join(results_dir, f"nhood_zscore_{sample_id}.tsv"), sep="\t"
        )

    # ── Co-occurrence ────────────────────────────────────────────────────
    sq.gr.co_occurrence(adata, cluster_key="cell_type")
    n_types = adata.obs["cell_type"].nunique()
    # ~2.5 inches per cell type, capped for readability
    fig_width = min(max(n_types * 2.5, 10), 60)
    fig_height = max(4, n_types * 0.6)

    sq.pl.co_occurrence(
        adata,
        cluster_key="cell_type",
        figsize=(fig_width, fig_height),
    )

    fig = plt.gcf()
    axes = fig.get_axes()

    seen_labels = set()
    for i, ax in enumerate(axes):
        # Remove duplicate legends — keep only the last non-empty one
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()

        # Rotate x-tick labels to prevent overlap
        ax.tick_params(axis="x", labelrotation=45)
        ax.tick_params(axis="y", labelsize=7)

        # Only show y-axis label on leftmost plots
        if i % n_types != 0:
            ax.set_ylabel("")

    # Collect legend handles from the last axis that has them
    handles, labels = [], []
    for ax in reversed(axes):
        h, l = ax.get_legend_handles_labels()
        if h:
            handles, labels = h, l
            break

    # Deduplicate legend entries
    unique = dict(zip(labels, handles))
    if unique:
        fig.legend(
            unique.values(),
            unique.keys(),
            loc="lower center",
            bbox_to_anchor=(0.5, -0.12),
            ncol=min(len(unique), 6),
            frameon=True,
            title="Cell type",
            fontsize=8,
            title_fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(
        os.path.join(results_dir, f"co_occurrence_{sample_id}.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Per-region analysis
    has_regions = (
        "region_annotation" in adata.obs.columns
        and adata.obs["region_annotation"].nunique() > 1
        and not all(adata.obs["region_annotation"] == "Unlabeled")
    )
    if has_regions:
        regions = [r for r in adata.obs["region_annotation"].unique()
                   if r not in ("Unlabeled", "Bubble")]
        for region in regions:
            rdata = adata[adata.obs["region_annotation"] == region].copy()
            if rdata.n_obs < 30 or rdata.obs["cell_type"].nunique() < 2:
                continue
            sq.gr.spatial_neighbors(rdata, coord_type="generic", delaunay=True)
            sq.gr.nhood_enrichment(rdata, cluster_key="cell_type")
            sq.pl.nhood_enrichment(rdata, cluster_key="cell_type", figsize=(8, 8))
            plt.title(f"Nhood enrichment – {region}")
            plt.savefig(os.path.join(results_dir, f"nhood_{region}_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()

        # ── Co-occurrence BY AREA (composite; separate from the global figure) ──
        # Per region, the self-co-occurrence of each cell type over distance (how
        # spatially clustered each type is within that area), one panel per region,
        # coloured by the shared cell-type palette.
        region_curves = {}
        for region in regions:
            rdata = adata[adata.obs["region_annotation"] == region].copy()
            rdata.obs["cell_type"] = rdata.obs["cell_type"].cat.remove_unused_categories()
            if rdata.n_obs < 30 or rdata.obs["cell_type"].nunique() < 2:
                continue
            try:
                sq.gr.co_occurrence(rdata, cluster_key="cell_type")
                occ = rdata.uns["cell_type_co_occurrence"]["occ"]
                interval = np.asarray(rdata.uns["cell_type_co_occurrence"]["interval"])
                mids = (interval[:-1] + interval[1:]) / 2.0
                cats = list(rdata.obs["cell_type"].cat.categories)
                region_curves[region] = (mids, {ct: occ[i, i, :] for i, ct in enumerate(cats)})
            except Exception as e:
                log.warning("  co-occurrence failed for region %s: %s", region, e)

        if region_curves:
            nreg = len(region_curves)
            ncol = min(3, nreg)
            nrow = int(np.ceil(nreg / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow), squeeze=False)
            for ax in axes.ravel():
                ax.set_axis_off()
            legend_handles = {}
            for k, (region, (mids, curves)) in enumerate(region_curves.items()):
                ax = axes[k // ncol][k % ncol]
                ax.set_axis_on()
                for ct, curve in curves.items():
                    color = cell_type_palette.get(str(ct), "#cccccc")
                    ax.plot(mids, curve, color=color, lw=1.3)
                    legend_handles.setdefault(
                        str(ct), plt.Line2D([0], [0], color=color, lw=1.3))
                ax.set_title(str(region), fontsize=9)
                ax.set_xlabel("distance")
                ax.set_ylabel("co-occurrence ratio")
            if legend_handles:
                fig.legend(list(legend_handles.values()), list(legend_handles.keys()),
                           loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
                           title="Cell type", fontsize=7, title_fontsize=8)
            fig.suptitle(f"Self co-occurrence by area — {sample_id}", fontsize=12, y=1.02)
            fig.tight_layout()
            fig.savefig(os.path.join(results_dir, f"co_occurrence_by_region_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close(fig)

    log.info("Neighbourhood analysis complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
