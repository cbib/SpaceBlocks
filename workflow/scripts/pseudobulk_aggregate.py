"""
pseudobulk_aggregate.py – Aggregate cells into pseudobulk matrices
====================================================================
Uses decoupler.get_pseudobulk for memory-efficient aggregation.
Saves count matrices and metadata as TSVs for downstream R-based DE.

Output structure per analysis level:
  matrices/
  ├── manifest.tsv              (lists all available subgroups)
  ├── counts_pooled.tsv         (by_region: all cell types pooled)
  ├── metadata_pooled.tsv
  ├── counts_{CellType}.tsv     (by_celltype_region: one per type)
  └── metadata_{CellType}.tsv
"""

import logging
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import decoupler as dc


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
log = logging.getLogger("pseudobulk_aggregate")

ANNOT_COL_MAP = {
    "tsv_annotation": "cell_type_tsv",
    "refined_annotation": "cell_type_refined",
    "ingest_annotation": "cell_type_ingest",
}


def _clean_nan(adata, col):
    mask = (
        adata.obs[col].notna()
        & (adata.obs[col].astype(str) != "nan")
        & (adata.obs[col].astype(str) != "")
    )
    n_drop = (~mask).sum()
    if n_drop > 0:
        log.info("  Filtering %d 'nan' entries from '%s'", n_drop, col)
    return adata[mask].copy()


def _save_pseudobulk(pdata, out_dir, prefix):
    """Save count matrix and metadata from a pseudobulk AnnData."""
    if sparse.issparse(pdata.X):
        counts = pd.DataFrame(pdata.X.toarray(), index=pdata.obs_names,
                              columns=pdata.var_names)
    else:
        counts = pd.DataFrame(np.asarray(pdata.X), index=pdata.obs_names,
                              columns=pdata.var_names)
    counts = counts.round().astype(int)
    counts.to_csv(os.path.join(out_dir, f"counts_{prefix}.tsv"), sep="\t")
    pdata.obs.to_csv(os.path.join(out_dir, f"metadata_{prefix}.tsv"), sep="\t")
    log.info("  Saved: %s (%d samples × %d genes)", prefix, pdata.n_obs, pdata.n_vars)


# ── Parameters ───────────────────────────────────────────────────────────────
annot_type     = snakemake.params.annot_type
analysis_level = snakemake.params.analysis_level
MIN_CELLS      = int(snakemake.params.min_cells_per_pseudobulk)
MIN_COUNTS     = int(snakemake.params.min_counts_per_pseudobulk)
agg_dir        = str(snakemake.output.agg_dir)

try:
    log.info("=" * 70)
    log.info("Pseudobulk aggregation: %s / %s", annot_type, analysis_level)
    log.info("  min_cells=%d, min_counts=%d", MIN_CELLS, MIN_COUNTS)
    log.info("=" * 70)

    os.makedirs(agg_dir, exist_ok=True)

    adata = sc.read_h5ad(str(snakemake.input.adata))
    log.info("Loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)

    cell_type_col = ANNOT_COL_MAP.get(annot_type)
    if cell_type_col and cell_type_col not in adata.obs.columns:
        log.warning("Column '%s' not found. Skipping.", cell_type_col)
        Path(os.path.join(agg_dir, f"SKIPPED_no_{cell_type_col}.txt")).write_text(
            f"Column {cell_type_col} not found.\n"
        )
        sys.exit(0)

    # Find sample column
    sample_col = None
    for c in ["sample", "sample_batch"]:
        if c in adata.obs.columns:
            sample_col = c
            break
    if sample_col is None:
        raise ValueError("No sample column found.")

    region_col = "region_annotation"
    if region_col not in adata.obs.columns:
        log.warning("No region_annotation. Skipping.")
        Path(os.path.join(agg_dir, "SKIPPED_no_regions.txt")).write_text("")
        sys.exit(0)

    # Clean NaN from grouping columns
    for col in [region_col, sample_col]:
        adata = _clean_nan(adata, col)
    if cell_type_col:
        adata = _clean_nan(adata, cell_type_col)

    # Filter regions
    valid_regions = [r for r in adata.obs[region_col].unique()
                     if r not in ("Unlabeled", "Bubble")]
    adata = adata[adata.obs[region_col].isin(valid_regions)].copy()
    log.info("Regions: %s (%d cells)", valid_regions, adata.n_obs)

    if len(valid_regions) < 2:
        Path(os.path.join(agg_dir, "SKIPPED_fewer_than_2_regions.txt")).write_text("")
        sys.exit(0)

    if "raw_counts" not in adata.layers:
        raise ValueError("raw_counts layer missing.")

    manifest_rows = []

    if analysis_level == "by_region":
        pdata = dc.get_pseudobulk(
            adata, sample_col=sample_col, groups_col=region_col,
            layer="raw_counts", mode="sum",
            min_cells=MIN_CELLS, min_counts=MIN_COUNTS,
        )
        _save_pseudobulk(pdata, agg_dir, "pooled")
        manifest_rows.append({
            "prefix": "pooled", "grouping": "all_cells",
            "n_samples": pdata.n_obs, "n_genes": pdata.n_vars,
            "condition_col": region_col, "sample_col": sample_col,
        })

    elif analysis_level == "by_celltype_region":
        adata = adata[adata.obs[cell_type_col] != "Unannotated"].copy()
        cell_types = sorted(adata.obs[cell_type_col].unique())

        for ct in cell_types:
            safe_ct = ct.replace("/", "_").replace(" ", "_")
            ct_adata = adata[adata.obs[cell_type_col] == ct].copy()

            if ct_adata.n_obs < MIN_CELLS * 2:
                log.warning("  '%s': too few cells (%d). Skipping.", ct, ct_adata.n_obs)
                continue

            try:
                pdata = dc.get_pseudobulk(
                    ct_adata, sample_col=sample_col, groups_col=region_col,
                    layer="raw_counts", mode="sum",
                    min_cells=MIN_CELLS, min_counts=MIN_COUNTS,
                )
            except Exception as e:
                log.warning("  '%s': aggregation failed: %s", ct, e)
                continue

            if pdata.n_obs < 4:
                log.warning("  '%s': too few pseudobulk samples (%d). Skipping.", ct, pdata.n_obs)
                continue

            _save_pseudobulk(pdata, agg_dir, safe_ct)
            manifest_rows.append({
                "prefix": safe_ct, "grouping": ct,
                "n_samples": pdata.n_obs, "n_genes": pdata.n_vars,
                "condition_col": region_col, "sample_col": sample_col,
            })

    elif analysis_level == "by_niche_region":
        log.warning("by_niche_region not yet implemented.")
        Path(os.path.join(agg_dir, "SKIPPED_niche_not_implemented.txt")).write_text("")

    # Save manifest
    if manifest_rows:
        pd.DataFrame(manifest_rows).to_csv(
            os.path.join(agg_dir, "manifest.tsv"), sep="\t", index=False
        )

    log.info("Aggregation complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
