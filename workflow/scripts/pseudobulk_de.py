"""
pseudobulk_de.py – Multi-level pseudobulk DE using decoupler + pyDESeq2
=========================================================================
Uses decoupler.get_pseudobulk for memory-efficient aggregation and
filter_by_expr for gene filtering (edgeR-style).

Improvements over previous version:
- Memory-efficient: decoupler aggregates without densifying full matrix
- Saves pseudobulk count matrices before DE for further exploration
- Filters 'nan' levels from all grouping categories
- Results ordered by abs(log2FoldChange) descending
- de_n_genes used for top-gene heatmaps per contrast
- decoupler volcano plots
"""

import itertools
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
from scipy import sparse
import decoupler as dc
from pydeseq2.dds import DeseqDataSet, DefaultInference
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


ANNOT_COL_MAP = {
    "tsv_annotation": "cell_type_tsv",
    "refined_annotation": "cell_type_refined",
    "ingest_annotation": "cell_type_ingest",
}


def _clean_nan_levels(adata, col):
    """Remove observations where col is NaN, 'nan', or empty string."""
    mask = (
        adata.obs[col].notna()
        & (adata.obs[col].astype(str) != "nan")
        & (adata.obs[col].astype(str) != "")
    )
    if (~mask).sum() > 0:
        log.info("    Filtering %d 'nan' entries from '%s'", (~mask).sum(), col)
    return adata[mask].copy()


def _save_pseudobulk_matrix(pdata, out_dir, prefix):
    """Save the pseudobulk count matrix and metadata as TSVs."""
    if sparse.issparse(pdata.X):
        counts = pd.DataFrame(pdata.X.toarray(), index=pdata.obs_names,
                              columns=pdata.var_names)
    else:
        counts = pd.DataFrame(pdata.X, index=pdata.obs_names,
                              columns=pdata.var_names)
    counts.to_csv(os.path.join(out_dir, f"{prefix}_counts.tsv"), sep="\t")
    pdata.obs.to_csv(os.path.join(out_dir, f"{prefix}_metadata.tsv"), sep="\t")
    log.info("    Saved pseudobulk matrix: %s (%d samples × %d genes)",
             prefix, pdata.n_obs, pdata.n_vars)


def _run_deseq2_contrast(pdata, condition_col, cond_a, cond_b, out_dir,
                         contrast_name, de_n_genes):
    """Run pyDESeq2 for one pairwise contrast, save results + plots."""
    log.info("      Contrast: %s", contrast_name)

    try:
        inference = DefaultInference(n_cpus=1)
        dds = DeseqDataSet(
            adata=pdata,
            design_factors=condition_col,
            ref_level=[condition_col, cond_b],
            refit_cooks=True,
            inference=inference,
        )
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=[condition_col, cond_a, cond_b],
                              inference=inference)
        stat_res.summary()
        results = stat_res.results_df.copy()

        # Order by absolute log2FoldChange (descending)
        results = results.reindex(
            results["log2FoldChange"].abs().sort_values(ascending=False).index
        )

        # Save full results
        results.to_csv(os.path.join(out_dir, f"{contrast_name}.tsv"), sep="\t")

        n_sig = (results["padj"] < 0.05).sum()
        n_up = ((results["padj"] < 0.05) & (results["log2FoldChange"] > 0)).sum()
        n_down = ((results["padj"] < 0.05) & (results["log2FoldChange"] < 0)).sum()
        log.info("      %d sig (↑%d, ↓%d)", n_sig, n_up, n_down)

        # Volcano plot using decoupler
        try:
            logFCs = results[["log2FoldChange"]].T.rename(
                index={"log2FoldChange": contrast_name}
            )
            pvals = results[["padj"]].T.rename(index={"padj": contrast_name})
            fig, ax = plt.subplots(figsize=(8, 6))
            dc.plot_volcano(logFCs, pvals, contrast_name, top=10,
                            sign_thr=0.05, lFCs_thr=0.5, ax=ax)
            plt.tight_layout()
            safe = contrast_name.replace(" ", "_")
            plt.savefig(os.path.join(out_dir, f"volcano_{safe}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            log.warning("      Volcano plot failed: %s", e)

        # Top-N DE gene heatmap (using pseudobulk expression)
        try:
            sig_genes = results[results["padj"] < 0.05]
            top_up = sig_genes[sig_genes["log2FoldChange"] > 0].head(de_n_genes).index.tolist()
            top_down = sig_genes[sig_genes["log2FoldChange"] < 0].head(de_n_genes).index.tolist()
            top_genes = top_up + top_down
            top_genes = [g for g in top_genes if g in pdata.var_names]

            if len(top_genes) >= 2:
                # Normalise pseudobulk for plotting
                pdata_norm = pdata.copy()
                sc.pp.normalize_total(pdata_norm, target_sum=1e4)
                sc.pp.log1p(pdata_norm)

                sc.pl.heatmap(
                    pdata_norm[:, top_genes],
                    var_names=top_genes,
                    groupby=condition_col,
                    swap_axes=True,
                    show_gene_labels=True,
                    figsize=(max(6, len(pdata_norm.obs_names) * 0.8), max(4, len(top_genes) * 0.4)),
                )
                plt.savefig(os.path.join(out_dir, f"heatmap_top{de_n_genes}_{safe}.png"),
                            dpi=300, bbox_inches="tight")
                plt.close()
        except Exception as e:
            log.warning("      Heatmap failed: %s", e)

        return {
            "contrast": contrast_name,
            "n_tested": len(results),
            "n_significant": n_sig,
            "n_up": n_up,
            "n_down": n_down,
        }

    except Exception as e:
        log.warning("      pyDESeq2 failed: %s", e)
        return {
            "contrast": contrast_name,
            "n_tested": 0,
            "n_significant": 0,
            "n_up": 0,
            "n_down": 0,
            "error": str(e),
        }


def _run_holistic_lrt(pdata, condition_col, out_dir):
    """LRT: full model (~ condition) vs reduced (~ 1)."""
    log.info("    Holistic LRT …")
    try:
        inference = DefaultInference(n_cpus=1)
        dds = DeseqDataSet(
            adata=pdata,
            design_factors=condition_col,
            refit_cooks=True,
            inference=inference,
        )
        dds.deseq2()
        stat_res = DeseqStats(dds, inference=inference)
        stat_res.summary()
        results = stat_res.results_df.copy()
        results = results.reindex(
            results["log2FoldChange"].abs().sort_values(ascending=False).index
        )
        results.to_csv(os.path.join(out_dir, "holistic_LRT.tsv"), sep="\t")
        n_sig = (results["padj"] < 0.05).sum()
        log.info("    LRT: %d significant genes (padj < 0.05)", n_sig)
        return results
    except Exception as e:
        log.warning("    LRT failed: %s", e)
        return None


# ── Parameters ───────────────────────────────────────────────────────────────
annot_type     = snakemake.params.annot_type
analysis_level = snakemake.params.analysis_level
MIN_CELLS      = int(snakemake.params.min_cells_per_pseudobulk)
MIN_REPS       = int(snakemake.params.min_replicates)
DE_N_GENES     = int(snakemake.params.de_n_genes)

adata_path     = str(snakemake.input.adata)
results_dir    = str(snakemake.output.results_dir)

try:
    log.info("=" * 70)
    log.info("Pseudobulk DE: annot_type=%s, level=%s", annot_type, analysis_level)
    log.info("  min_cells=%d, min_replicates=%d, de_n_genes=%d",
             MIN_CELLS, MIN_REPS, DE_N_GENES)
    log.info("=" * 70)

    os.makedirs(results_dir, exist_ok=True)

    adata = sc.read_h5ad(adata_path)
    log.info("Loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # Determine cell-type column
    cell_type_col = ANNOT_COL_MAP.get(annot_type)
    if cell_type_col and cell_type_col not in adata.obs.columns:
        log.warning("Column '%s' not found. Skipping.", cell_type_col)
        Path(os.path.join(results_dir, f"SKIPPED_no_{cell_type_col}.txt")).write_text(
            f"Column {cell_type_col} not found.\n"
        )
        sys.exit(0)

    # Determine sample column
    sample_col = None
    for candidate in ["sample", "sample_batch"]:
        if candidate in adata.obs.columns:
            sample_col = candidate
            break
    if sample_col is None:
        raise ValueError("No sample column found in adata.obs")

    # Check region annotations
    region_col = "region_annotation"
    if region_col not in adata.obs.columns:
        log.warning("No region_annotation. Skipping.")
        Path(os.path.join(results_dir, "SKIPPED_no_regions.txt")).write_text(
            "No region_annotation column.\n"
        )
        sys.exit(0)

    # Clean nan levels from all grouping columns
    for col in [region_col, sample_col]:
        adata = _clean_nan_levels(adata, col)
    if cell_type_col:
        adata = _clean_nan_levels(adata, cell_type_col)

    # Filter regions
    valid_regions = [r for r in adata.obs[region_col].unique()
                     if r not in ("Unlabeled", "Bubble")]
    if len(valid_regions) < 2:
        log.warning("Fewer than 2 valid regions. Skipping.")
        Path(os.path.join(results_dir, "SKIPPED_fewer_than_2_regions.txt")).write_text(
            f"Only {len(valid_regions)} region(s).\n"
        )
        sys.exit(0)

    adata = adata[adata.obs[region_col].isin(valid_regions)].copy()
    log.info("After filtering: %d cells, %d regions: %s",
             adata.n_obs, len(valid_regions), valid_regions)

    # Ensure raw_counts layer exists
    if "raw_counts" not in adata.layers:
        raise ValueError("raw_counts layer missing.")

    # ── Dispatch by analysis level ───────────────────────────────────────
    if analysis_level == "by_region":
        log.info("Level: by_region (all cell types pooled)")

        # Pseudobulk: one profile per sample × region
        pdata = dc.get_pseudobulk(
            adata,
            sample_col=sample_col,
            groups_col=region_col,
            layer="raw_counts",
            mode="sum",
            min_cells=MIN_CELLS,
            min_counts=1000,
        )
        log.info("  Pseudobulk: %d samples × %d genes", pdata.n_obs, pdata.n_vars)

        # Save matrix
        _save_pseudobulk_matrix(pdata, results_dir, "pseudobulk_by_region")

        # Filter genes (edgeR-style)
        dc.filter_by_expr(pdata, group=region_col, min_count=10, min_total_count=15)
        log.info("  After gene filter: %d genes", pdata.n_vars)

        # Check replicates
        region_counts = pdata.obs[region_col].value_counts()
        valid_conds = region_counts[region_counts >= MIN_REPS].index
        if len(valid_conds) < 2:
            log.warning("  < 2 regions with >= %d replicates. Skipping.", MIN_REPS)
            Path(os.path.join(results_dir, "SKIPPED_insufficient.txt")).write_text(
                "Insufficient replicates.\n"
            )
            sys.exit(0)

        pdata = pdata[pdata.obs[region_col].isin(valid_conds)].copy()

        # LRT
        _run_holistic_lrt(pdata, region_col, results_dir)

        # Pairwise
        pairwise_dir = os.path.join(results_dir, "pairwise")
        os.makedirs(pairwise_dir, exist_ok=True)
        conditions = sorted(valid_conds)
        summary_rows = []
        for cond_a, cond_b in itertools.combinations(conditions, 2):
            contrast_name = f"{cond_a}_vs_{cond_b}"
            sub = pdata[pdata.obs[region_col].isin([cond_a, cond_b])].copy()
            row = _run_deseq2_contrast(
                sub, region_col, cond_a, cond_b,
                pairwise_dir, contrast_name, DE_N_GENES,
            )
            summary_rows.append(row)
        if summary_rows:
            pd.DataFrame(summary_rows).to_csv(
                os.path.join(results_dir, "pairwise_summary.tsv"),
                sep="\t", index=False,
            )

    elif analysis_level == "by_celltype_region":
        log.info("Level: by_celltype_region")

        if cell_type_col is None:
            log.warning("No cell type column. Skipping.")
            sys.exit(0)

        # Filter unannotated
        adata = adata[adata.obs[cell_type_col] != "Unannotated"].copy()
        cell_types = sorted(adata.obs[cell_type_col].unique())
        log.info("  Cell types: %s", cell_types)

        for ct in cell_types:
            safe_ct = ct.replace("/", "_").replace(" ", "_")
            ct_dir = os.path.join(results_dir, safe_ct)
            os.makedirs(ct_dir, exist_ok=True)

            ct_adata = adata[adata.obs[cell_type_col] == ct].copy()
            log.info("  Cell type '%s': %d cells", ct, ct_adata.n_obs)

            if ct_adata.n_obs < MIN_CELLS * 2:
                log.warning("    Too few cells. Skipping.")
                Path(os.path.join(ct_dir, "SKIPPED_too_few_cells.txt")).write_text(
                    f"{ct_adata.n_obs} cells.\n"
                )
                continue

            # Pseudobulk per sample × region for this cell type
            try:
                pdata = dc.get_pseudobulk(
                    ct_adata,
                    sample_col=sample_col,
                    groups_col=region_col,
                    layer="raw_counts",
                    mode="sum",
                    min_cells=MIN_CELLS,
                    min_counts=1000,
                )
            except Exception as e:
                log.warning("    Pseudobulk aggregation failed: %s", e)
                continue

            if pdata.n_obs < 4:
                log.warning("    Too few pseudobulk samples (%d). Skipping.", pdata.n_obs)
                continue

            # Save matrix
            _save_pseudobulk_matrix(pdata, ct_dir, f"pseudobulk_{safe_ct}")

            # Filter genes
            try:
                dc.filter_by_expr(pdata, group=region_col, min_count=10,
                                  min_total_count=15)
            except Exception as e:
                log.warning("    Gene filtering failed: %s", e)

            # Check replicates
            region_counts = pdata.obs[region_col].value_counts()
            valid_conds = region_counts[region_counts >= MIN_REPS].index
            if len(valid_conds) < 2:
                log.warning("    < 2 regions with >= %d replicates. Skipping.", MIN_REPS)
                Path(os.path.join(ct_dir, "SKIPPED_insufficient.txt")).write_text(
                    "Insufficient replicates.\n"
                )
                continue

            pdata = pdata[pdata.obs[region_col].isin(valid_conds)].copy()
            log.info("    %d pseudobulk samples, %d regions, %d genes",
                     pdata.n_obs, len(valid_conds), pdata.n_vars)

            # LRT
            _run_holistic_lrt(pdata, region_col, ct_dir)

            # Pairwise
            pairwise_dir = os.path.join(ct_dir, "pairwise")
            os.makedirs(pairwise_dir, exist_ok=True)
            conditions = sorted(valid_conds)
            summary_rows = []
            for cond_a, cond_b in itertools.combinations(conditions, 2):
                contrast_name = f"{cond_a}_vs_{cond_b}"
                sub = pdata[pdata.obs[region_col].isin([cond_a, cond_b])].copy()
                row = _run_deseq2_contrast(
                    sub, region_col, cond_a, cond_b,
                    pairwise_dir, contrast_name, DE_N_GENES,
                )
                summary_rows.append(row)
            if summary_rows:
                pd.DataFrame(summary_rows).to_csv(
                    os.path.join(ct_dir, "pairwise_summary.tsv"),
                    sep="\t", index=False,
                )

    elif analysis_level == "by_niche_region":
        log.warning("by_niche_region not yet implemented.")
        Path(os.path.join(results_dir, "SKIPPED_niche_not_implemented.txt")).write_text(
            "Niche identification must run first.\n"
        )

    else:
        raise ValueError(f"Unknown analysis_level: {analysis_level}")

    log.info("Pseudobulk DE complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
