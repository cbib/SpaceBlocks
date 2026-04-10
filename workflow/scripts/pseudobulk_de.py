"""
pseudobulk_de.py – Pseudobulk DE with pyDESeq2
=================================================
Reads the annotated adata and uses the cell_type column specified
by the annot_type param ('tsv_annotation' → 'cell_type_tsv',
'refined_annotation' → 'cell_type_refined').
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

sample_id   = snakemake.params.sample_id
annot_type  = snakemake.params.annot_type
adata_path  = str(snakemake.input.adata)
results_dir = str(snakemake.output.results_dir)

# Map annot_type to obs column
ANNOT_COL_MAP = {
    "tsv_annotation": "cell_type_tsv",
    "refined_annotation": "cell_type_refined",
}

try:
    log.info("=" * 70)
    log.info("Pseudobulk DE: sample=%s, annot_type=%s", sample_id, annot_type)
    log.info("=" * 70)

    os.makedirs(results_dir, exist_ok=True)

    cell_type_col = ANNOT_COL_MAP.get(annot_type)
    if cell_type_col is None:
        raise ValueError(f"Unknown annot_type '{annot_type}'. Expected: {list(ANNOT_COL_MAP.keys())}")

    adata = sc.read_h5ad(adata_path)

    if cell_type_col not in adata.obs.columns:
        log.warning("Column '%s' not found in adata. Skipping.", cell_type_col)
        Path(os.path.join(results_dir, f"SKIPPED_no_{cell_type_col}.txt")).write_text(
            f"Column {cell_type_col} not found.\n"
        )
        sys.exit(0)

    log.info("Using annotation column: %s", cell_type_col)
    adata.obs["cell_type"] = adata.obs[cell_type_col]

    # Check regions
    has_regions = (
        "region_annotation" in adata.obs.columns
        and adata.obs["region_annotation"].nunique() > 1
        and not all(adata.obs["region_annotation"] == "Unlabeled")
    )

    if not has_regions:
        log.warning("No region annotations. Skipping DE.")
        Path(os.path.join(results_dir, "SKIPPED_no_regions.txt")).write_text(
            "No region annotations found.\n"
        )
        sys.exit(0)

    regions = [r for r in adata.obs["region_annotation"].unique()
               if r not in ("Unlabeled", "Bubble")]
    if len(regions) < 2:
        log.warning("Fewer than 2 regions. Skipping.")
        Path(os.path.join(results_dir, "SKIPPED_fewer_than_2_regions.txt")).write_text(
            f"Only {len(regions)} region(s).\n"
        )
        sys.exit(0)

    mask = (
        (adata.obs["cell_type"] != "Unannotated")
        & adata.obs["region_annotation"].isin(regions)
    )
    adata_sub = adata[mask].copy()

    if "raw_counts" not in adata_sub.layers:
        log.error("raw_counts layer missing.")
        raise ValueError("raw_counts layer required for pseudobulk DE.")

    for ct in adata_sub.obs["cell_type"].unique():
        if ct == "Unannotated":
            continue
        log.info("DE for cell type: %s", ct)
        ct_adata = adata_sub[adata_sub.obs["cell_type"] == ct].copy()

        if ct_adata.n_obs < 10:
            log.warning("  Too few cells (%d). Skipping.", ct_adata.n_obs)
            continue

        aggregated = sc.get.aggregate(ct_adata, by="region_annotation", func=["sum"], layer="raw_counts")
        if aggregated.n_obs < 2:
            continue

        counts_df = pd.DataFrame(
            aggregated.layers["sum"], index=aggregated.obs_names, columns=aggregated.var_names
        ).astype(int)
        metadata_df = pd.DataFrame(
            {"region": aggregated.obs["region_annotation"].values}, index=aggregated.obs_names
        )
        counts_df = counts_df.loc[:, counts_df.sum(axis=0) > 0]
        if counts_df.shape[1] < 10:
            continue

        try:
            dds = DeseqDataSet(counts=counts_df, metadata=metadata_df, design="~region", refit_cooks=True)
            dds.deseq2()
            region_vals = metadata_df["region"].unique()
            for i, ref in enumerate(region_vals):
                for test in region_vals[i + 1:]:
                    ds = DeseqStats(dds, contrast=["region", test, ref])
                    ds.summary()
                    out_file = os.path.join(results_dir, f"DE_{ct}_{test}_vs_{ref}.tsv")
                    ds.results_df.to_csv(out_file, sep="\t")
                    log.info("  Saved: %s", out_file)
        except Exception as e:
            log.warning("  pyDESeq2 failed for '%s': %s", ct, e)

    log.info("Pseudobulk DE complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
