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
adata_path  = str(snakemake.input.adata)
results_dir = str(snakemake.output.results_dir)

ANNOT_COL_MAP = {
    "tsv_annotation": "cell_type_tsv",
    "refined_annotation": "cell_type_refined",
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

    log.info("Neighbourhood analysis complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
