"""
pseudobulk_de.py – Pseudobulk differential expression with pyDESeq2
=====================================================================
Reads the user-provided cluster annotation TSV, maps Leiden clusters to
cell-type labels, then performs pseudobulk aggregation and differential
expression per cell type across region annotations using pyDESeq2.

The annotation TSV format:
- Rows: cluster numbers (0, 1, 2, …)
- Columns: sample names
- Values: cell-type labels assigned by the user

For each cell type that has cells in at least 2 regions, a DESeq2
comparison is run across regions.  Results are saved as TSV files.
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
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


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
log = logging.getLogger("pseudobulk_de")

# ── Parameters ───────────────────────────────────────────────────────────────
sample_id         = snakemake.params.sample_id
adata_path        = str(snakemake.input.adata)
annotations_path  = str(snakemake.input.cluster_annotations)
results_dir       = str(snakemake.output.results_dir)

try:
    log.info("=" * 70)
    log.info("Pseudobulk DE for sample: %s", sample_id)
    log.info("=" * 70)

    os.makedirs(results_dir, exist_ok=True)

    # ── Load annotation TSV and build cluster→cell_type mapping ──────────
    log.info("Loading cluster annotations from: %s", annotations_path)
    annot_df = pd.read_csv(annotations_path, sep="\t", index_col=0)

    if sample_id not in annot_df.columns:
        log.warning("Sample '%s' not found in annotation TSV columns: %s. Skipping.",
                     sample_id, list(annot_df.columns))
        Path(os.path.join(results_dir, "SKIPPED_sample_not_in_annotations.txt")).write_text(
            f"Sample {sample_id} not found in {annotations_path}\n"
        )
        sys.exit(0)

    cluster_to_celltype = annot_df[sample_id].dropna().to_dict()
    cluster_to_celltype = {str(k): v for k, v in cluster_to_celltype.items()}
    log.info("Cluster → cell type mapping: %s", cluster_to_celltype)

    # ── Load adata ───────────────────────────────────────────────────────
    log.info("Loading adata: %s", adata_path)
    adata = sc.read_h5ad(adata_path)

    if "leiden" not in adata.obs.columns:
        raise ValueError("adata.obs missing 'leiden' column — run post_processing first.")
    if "raw_counts" not in adata.layers:
        raise ValueError("adata.layers missing 'raw_counts' — ensure post_processing saved it.")

    # ── Map clusters to cell types ───────────────────────────────────────
    adata.obs["cell_type"] = (
        adata.obs["leiden"].astype(str).map(cluster_to_celltype).fillna("Unannotated")
    )
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    log.info("Cell type distribution:\n%s", adata.obs["cell_type"].value_counts().to_string())

    # Save cell type annotation plot
    if "leiden" in adata.obs.columns:
        ct_counts = pd.crosstab(adata.obs["leiden"], adata.obs["cell_type"])
        ct_counts.plot(kind="bar", stacked=True, figsize=(12, 6))
        plt.title(f"Cell type composition by cluster – {sample_id}")
        plt.ylabel("Number of cells")
        plt.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.savefig(os.path.join(results_dir, f"celltype_by_cluster_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    # ── Check if region annotations are available ────────────────────────
    has_regions = (
        "region_annotation" in adata.obs.columns
        and adata.obs["region_annotation"].nunique() > 1
        and not all(adata.obs["region_annotation"].str.contains("not_found", case=False, na=False))
    )

    if not has_regions:
        log.warning("No region annotations available. Skipping pseudobulk DE.")
        Path(os.path.join(results_dir, "SKIPPED_no_region_annotations.txt")).write_text(
            "No region annotations found in adata.obs['region_annotation'].\n"
            "Pseudobulk DE requires at least 2 regions to compare.\n"
        )
        sys.exit(0)

    regions = adata.obs["region_annotation"].unique()
    regions = [r for r in regions if r not in ("Unlabeled", "Bubble")]
    log.info("Regions for DE: %s", regions)

    if len(regions) < 2:
        log.warning("Fewer than 2 annotated regions. Skipping DE.")
        Path(os.path.join(results_dir, "SKIPPED_fewer_than_2_regions.txt")).write_text(
            f"Only {len(regions)} annotated regions found. Need at least 2.\n"
        )
        sys.exit(0)

    # ── Pseudobulk DE per cell type ──────────────────────────────────────
    # Filter to annotated cells and regions
    mask = (
        (adata.obs["cell_type"] != "Unannotated")
        & (adata.obs["region_annotation"].isin(regions))
    )
    adata_sub = adata[mask].copy()

    cell_types = [ct for ct in adata_sub.obs["cell_type"].unique() if ct != "Unannotated"]

    for ct in cell_types:
        log.info("Running pseudobulk DE for cell type: %s", ct)

        ct_adata = adata_sub[adata_sub.obs["cell_type"] == ct].copy()

        if ct_adata.n_obs < 10:
            log.warning("  Too few cells (%d) for cell type '%s'. Skipping.", ct_adata.n_obs, ct)
            continue

        # Aggregate raw counts by region (pseudobulk)
        aggregated = sc.get.aggregate(
            ct_adata, by="region_annotation", func=["sum"], layer="raw_counts"
        )

        if aggregated.n_obs < 2:
            log.warning("  Fewer than 2 pseudobulk samples for '%s'. Skipping.", ct)
            continue

        # Build counts and metadata for pyDESeq2
        counts_df = pd.DataFrame(
            aggregated.layers["sum"],
            index=aggregated.obs_names,
            columns=aggregated.var_names,
        ).astype(int)

        metadata_df = pd.DataFrame(
            {"region": aggregated.obs["region_annotation"].values},
            index=aggregated.obs_names,
        )

        # Filter genes with zero total counts
        counts_df = counts_df.loc[:, counts_df.sum(axis=0) > 0]

        if counts_df.shape[1] < 10:
            log.warning("  Too few expressed genes (%d) for '%s'. Skipping.", counts_df.shape[1], ct)
            continue

        try:
            dds = DeseqDataSet(
                counts=counts_df,
                metadata=metadata_df,
                design="~region",
                refit_cooks=True,
            )
            dds.deseq2()

            # Compare each pair of regions (first vs rest for simplicity)
            region_vals = metadata_df["region"].unique()
            for i, ref_region in enumerate(region_vals):
                for test_region in region_vals[i + 1:]:
                    log.info("  Contrast: %s vs %s", test_region, ref_region)
                    ds = DeseqStats(dds, contrast=["region", test_region, ref_region])
                    ds.summary()

                    out_file = os.path.join(
                        results_dir,
                        f"DE_{ct}_{test_region}_vs_{ref_region}.tsv"
                    )
                    ds.results_df.to_csv(out_file, sep="\t")
                    log.info("  Saved: %s", out_file)

        except Exception as e:
            log.warning("  pyDESeq2 failed for '%s': %s", ct, e)
            continue

    # Save the annotated adata with cell types
    adata.write(os.path.join(results_dir, f"adata_annotated_{sample_id}.h5ad"))

    log.info("Pseudobulk DE complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
