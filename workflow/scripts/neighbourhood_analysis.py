"""
neighbourhood_analysis.py – Spatial neighbourhood analysis with squidpy
=========================================================================
Computes neighbourhood enrichment and co-occurrence scores based on
user-provided cell-type annotations.  Uses squidpy's spatial graph
functions on the cell centroid coordinates stored in adata.obsm["spatial"].
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
log = logging.getLogger("neighbourhood_analysis")

# ── Parameters ───────────────────────────────────────────────────────────────
sample_id         = snakemake.params.sample_id
adata_path        = str(snakemake.input.adata)
annotations_path  = str(snakemake.input.cluster_annotations)
results_dir       = str(snakemake.output.results_dir)

try:
    log.info("=" * 70)
    log.info("Neighbourhood analysis for sample: %s", sample_id)
    log.info("=" * 70)

    os.makedirs(results_dir, exist_ok=True)

    # ── Load annotation TSV ──────────────────────────────────────────────
    annot_df = pd.read_csv(annotations_path, sep="\t", index_col=0)

    if sample_id not in annot_df.columns:
        log.warning("Sample '%s' not found in annotation TSV. Skipping.", sample_id)
        Path(os.path.join(results_dir, "SKIPPED_sample_not_in_annotations.txt")).write_text(
            f"Sample {sample_id} not found in {annotations_path}\n"
        )
        sys.exit(0)

    cluster_to_celltype = {str(k): v for k, v in annot_df[sample_id].dropna().to_dict().items()}

    # ── Load adata ───────────────────────────────────────────────────────
    log.info("Loading adata: %s", adata_path)
    adata = sc.read_h5ad(adata_path)

    # Map clusters to cell types
    adata.obs["cell_type"] = (
        adata.obs["leiden"].astype(str).map(cluster_to_celltype).fillna("Unannotated")
    )
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    # Filter out unannotated cells
    adata = adata[adata.obs["cell_type"] != "Unannotated"].copy()

    if adata.n_obs < 50:
        log.warning("Too few annotated cells (%d). Skipping.", adata.n_obs)
        Path(os.path.join(results_dir, "SKIPPED_too_few_cells.txt")).write_text(
            f"Only {adata.n_obs} annotated cells found.\n"
        )
        sys.exit(0)

    n_cell_types = adata.obs["cell_type"].nunique()
    if n_cell_types < 2:
        log.warning("Fewer than 2 cell types (%d). Skipping.", n_cell_types)
        Path(os.path.join(results_dir, "SKIPPED_fewer_than_2_celltypes.txt")).write_text(
            f"Only {n_cell_types} cell type(s) found.\n"
        )
        sys.exit(0)

    log.info("Annotated cells: %d, cell types: %d", adata.n_obs, n_cell_types)

    # ── Build spatial graph ──────────────────────────────────────────────
    log.info("Building spatial neighbours graph (Delaunay) …")
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    # ── Neighbourhood enrichment ─────────────────────────────────────────
    log.info("Computing neighbourhood enrichment …")
    sq.gr.nhood_enrichment(adata, cluster_key="cell_type")

    sq.pl.nhood_enrichment(adata, cluster_key="cell_type", figsize=(8, 8))
    plt.savefig(os.path.join(results_dir, f"nhood_enrichment_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # Save enrichment z-scores
    if "cell_type_nhood_enrichment" in adata.uns:
        zscore = adata.uns["cell_type_nhood_enrichment"]["zscore"]
        categories = adata.obs["cell_type"].cat.categories
        zscore_df = pd.DataFrame(zscore, index=categories, columns=categories)
        zscore_df.to_csv(os.path.join(results_dir, f"nhood_enrichment_zscore_{sample_id}.tsv"),
                         sep="\t")

    # ── Co-occurrence ────────────────────────────────────────────────────
    log.info("Computing co-occurrence …")
    sq.gr.co_occurrence(adata, cluster_key="cell_type")

    sq.pl.co_occurrence(adata, cluster_key="cell_type", figsize=(10, 6))
    plt.savefig(os.path.join(results_dir, f"co_occurrence_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # ── Per-region neighbourhood analysis (if region annotations exist) ──
    has_regions = (
        "region_annotation" in adata.obs.columns
        and adata.obs["region_annotation"].nunique() > 1
        and not all(adata.obs["region_annotation"].str.contains("not_found", case=False, na=False))
    )

    if has_regions:
        regions = [r for r in adata.obs["region_annotation"].unique()
                   if r not in ("Unlabeled", "Bubble")]
        log.info("Running per-region neighbourhood analysis for: %s", regions)

        for region in regions:
            region_adata = adata[adata.obs["region_annotation"] == region].copy()
            if region_adata.n_obs < 30 or region_adata.obs["cell_type"].nunique() < 2:
                log.warning("  Region '%s': too few cells/types. Skipping.", region)
                continue

            sq.gr.spatial_neighbors(region_adata, coord_type="generic", delaunay=True)
            sq.gr.nhood_enrichment(region_adata, cluster_key="cell_type")

            sq.pl.nhood_enrichment(region_adata, cluster_key="cell_type", figsize=(8, 8))
            plt.title(f"Neighbourhood enrichment – {region}")
            plt.savefig(
                os.path.join(results_dir, f"nhood_enrichment_{region}_{sample_id}.png"),
                dpi=300, bbox_inches="tight",
            )
            plt.close()

            if "cell_type_nhood_enrichment" in region_adata.uns:
                zscore = region_adata.uns["cell_type_nhood_enrichment"]["zscore"]
                categories = region_adata.obs["cell_type"].cat.categories
                zscore_df = pd.DataFrame(zscore, index=categories, columns=categories)
                zscore_df.to_csv(
                    os.path.join(results_dir, f"nhood_enrichment_zscore_{region}_{sample_id}.tsv"),
                    sep="\t",
                )

    log.info("Neighbourhood analysis complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
