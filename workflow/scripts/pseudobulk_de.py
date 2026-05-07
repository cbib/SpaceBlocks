"""
pseudobulk_de.py – Multi-level pseudobulk differential expression
====================================================================
Aggregates single cells into pseudobulk profiles (one per sample × group)
and runs pyDESeq2 for DE analysis.

Analysis levels:
- by_region:           aggregate per sample × region, test region effect
- by_celltype_region:  aggregate per sample × region × cell_type,
                       test region effect within each cell type
- by_niche_region:     (future) same, using niche labels

For each level, two test types:
- Holistic LRT:  does the factor explain variance? (full vs reduced model)
- Pairwise Wald: all unique region pairs, BH-corrected
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
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from statsmodels.stats.multitest import multipletests


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


# ── Helpers ──────────────────────────────────────────────────────────────────

ANNOT_COL_MAP = {
    "tsv_annotation": "cell_type_tsv",
    "refined_annotation": "cell_type_refined",
    "ingest_annotation": "cell_type_ingest",
}


def get_raw_counts(adata):
    """Extract raw integer counts, preferring the raw_counts layer."""
    if "raw_counts" in adata.layers:
        X = adata.layers["raw_counts"]
    else:
        X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    return np.round(X).astype(int)


def aggregate_pseudobulk(adata, groupby_cols):
    """
    Sum raw counts per group defined by groupby_cols.
    Returns (counts_df, metadata_df) where each row is one pseudobulk sample.
    """
    raw = get_raw_counts(adata)

    # Build group keys
    group_key = adata.obs[groupby_cols].astype(str).agg("__".join, axis=1)
    unique_groups = group_key.unique()

    counts_list = []
    meta_list = []

    for grp in unique_groups:
        mask = group_key == grp
        n_cells = mask.sum()
        summed = raw[mask.values].sum(axis=0)
        counts_list.append(summed)
        parts = dict(zip(groupby_cols, grp.split("__")))
        parts["n_cells"] = n_cells
        parts["pseudobulk_id"] = grp
        meta_list.append(parts)

    counts_df = pd.DataFrame(
        np.array(counts_list),
        columns=adata.var_names,
        index=[m["pseudobulk_id"] for m in meta_list],
    )
    metadata_df = pd.DataFrame(meta_list).set_index("pseudobulk_id")

    return counts_df, metadata_df


def filter_pseudobulk(counts_df, metadata_df, min_cells, min_replicates,
                      condition_col):
    """
    Filter pseudobulk samples:
    - Remove samples with fewer than min_cells
    - Remove conditions with fewer than min_replicates samples
    - Remove genes with zero total counts
    """
    # Min cells
    keep = metadata_df["n_cells"] >= min_cells
    if (~keep).sum() > 0:
        log.info("    Dropping %d pseudobulk samples with < %d cells",
                 (~keep).sum(), min_cells)
    counts_df = counts_df[keep]
    metadata_df = metadata_df[keep]

    # Min replicates per condition
    cond_counts = metadata_df[condition_col].value_counts()
    valid_conds = cond_counts[cond_counts >= min_replicates].index
    keep = metadata_df[condition_col].isin(valid_conds)
    if (~keep).sum() > 0:
        dropped_conds = set(metadata_df[condition_col]) - set(valid_conds)
        log.info("    Dropping conditions with < %d replicates: %s",
                 min_replicates, dropped_conds)
    counts_df = counts_df[keep]
    metadata_df = metadata_df[keep]

    # Remove zero-count genes
    nonzero = counts_df.sum(axis=0) > 0
    counts_df = counts_df.loc[:, nonzero]

    return counts_df, metadata_df


def run_holistic_lrt(counts_df, metadata_df, condition_col, out_dir):
    """
    Likelihood ratio test: full model (~ condition) vs reduced (~ 1).
    Tests whether condition has ANY effect on gene expression.
    """
    log.info("    Running holistic LRT …")
    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design=f"~ {condition_col}",
            refit_cooks=True,
        )
        dds.deseq2()

        # LRT
        stat_res = DeseqStats(dds, contrast=None)
        stat_res.summary()
        results = stat_res.results_df.copy()
        results = results.sort_values("padj")
        results.to_csv(os.path.join(out_dir, "holistic_LRT.tsv"), sep="\t")

        n_sig = (results["padj"] < 0.05).sum()
        log.info("    LRT: %d significant genes (padj < 0.05)", n_sig)

        return results
    except Exception as e:
        log.warning("    LRT failed: %s", e)
        return None


def run_pairwise_wald(counts_df, metadata_df, condition_col, out_dir,
                      de_n_genes):
    """
    Run Wald test for all pairwise contrasts of the condition factor.
    Each contrast is BH-corrected independently. A summary across all
    contrasts is also saved.
    """
    pairwise_dir = os.path.join(out_dir, "pairwise")
    os.makedirs(pairwise_dir, exist_ok=True)

    conditions = sorted(metadata_df[condition_col].unique())
    pairs = list(itertools.combinations(conditions, 2))
    log.info("    Running %d pairwise Wald tests: %s", len(pairs), pairs)

    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design=f"~ {condition_col}",
            refit_cooks=True,
        )
        dds.deseq2()
    except Exception as e:
        log.warning("    pyDESeq2 fitting failed: %s", e)
        return

    summary_rows = []

    for cond_a, cond_b in pairs:
        contrast_name = f"{cond_a}_vs_{cond_b}"
        log.info("      Contrast: %s", contrast_name)
        try:
            stat_res = DeseqStats(
                dds, contrast=[condition_col, cond_a, cond_b]
            )
            stat_res.summary()
            results = stat_res.results_df.copy()
            results = results.sort_values("padj")

            # Save full results
            results.to_csv(
                os.path.join(pairwise_dir, f"{contrast_name}.tsv"), sep="\t"
            )

            n_sig = (results["padj"] < 0.05).sum()
            n_up = ((results["padj"] < 0.05) & (results["log2FoldChange"] > 0)).sum()
            n_down = ((results["padj"] < 0.05) & (results["log2FoldChange"] < 0)).sum()
            summary_rows.append({
                "contrast": contrast_name,
                "n_tested": len(results),
                "n_significant": n_sig,
                "n_up": n_up,
                "n_down": n_down,
            })
            log.info("      %d sig (↑%d, ↓%d)", n_sig, n_up, n_down)

            # Volcano plot
            _volcano_plot(results, contrast_name, pairwise_dir)

            # MA plot
            _ma_plot(results, contrast_name, pairwise_dir)

        except Exception as e:
            log.warning("      Contrast %s failed: %s", contrast_name, e)
            summary_rows.append({
                "contrast": contrast_name,
                "n_tested": 0,
                "n_significant": 0,
                "n_up": 0,
                "n_down": 0,
                "error": str(e),
            })

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            os.path.join(out_dir, "pairwise_summary.tsv"), sep="\t", index=False
        )


def _volcano_plot(results, title, out_dir):
    """Standard volcano plot: log2FC vs -log10(padj)."""
    fig, ax = plt.subplots(figsize=(8, 6))
    r = results.dropna(subset=["padj", "log2FoldChange"])
    neg_log_p = -np.log10(r["padj"].clip(lower=1e-300))

    sig = r["padj"] < 0.05
    up = sig & (r["log2FoldChange"] > 1)
    down = sig & (r["log2FoldChange"] < -1)
    ns = ~(up | down)

    ax.scatter(r.loc[ns, "log2FoldChange"], neg_log_p[ns],
               c="#aaaaaa", s=3, alpha=0.5, label="NS")
    ax.scatter(r.loc[up, "log2FoldChange"], neg_log_p[up],
               c="#e41a1c", s=5, alpha=0.7, label=f"Up ({up.sum()})")
    ax.scatter(r.loc[down, "log2FoldChange"], neg_log_p[down],
               c="#377eb8", s=5, alpha=0.7, label=f"Down ({down.sum()})")

    ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=0.5)
    ax.axvline(1, ls="--", c="grey", lw=0.5)
    ax.axvline(-1, ls="--", c="grey", lw=0.5)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    plt.tight_layout()
    safe = title.replace(" ", "_")
    plt.savefig(os.path.join(out_dir, f"volcano_{safe}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()


def _ma_plot(results, title, out_dir):
    """MA plot: log2FC vs mean expression (baseMean)."""
    fig, ax = plt.subplots(figsize=(8, 6))
    r = results.dropna(subset=["padj", "log2FoldChange", "baseMean"])
    r = r[r["baseMean"] > 0]

    sig = r["padj"] < 0.05
    ax.scatter(np.log10(r.loc[~sig, "baseMean"]), r.loc[~sig, "log2FoldChange"],
               c="#aaaaaa", s=3, alpha=0.5, label="NS")
    ax.scatter(np.log10(r.loc[sig, "baseMean"]), r.loc[sig, "log2FoldChange"],
               c="#e41a1c", s=5, alpha=0.7, label=f"Sig ({sig.sum()})")

    ax.axhline(0, ls="-", c="black", lw=0.5)
    ax.set_xlabel("log10(baseMean)")
    ax.set_ylabel("log2 Fold Change")
    ax.set_title(title)
    ax.legend(fontsize=8)
    plt.tight_layout()
    safe = title.replace(" ", "_")
    plt.savefig(os.path.join(out_dir, f"MA_{safe}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()


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
    log.info("  min_cells=%d, min_replicates=%d", MIN_CELLS, MIN_REPS)
    log.info("=" * 70)

    os.makedirs(results_dir, exist_ok=True)

    adata = sc.read_h5ad(adata_path)
    log.info("Loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # Determine which cell-type column to use
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
        log.warning("No region_annotation column. Skipping.")
        Path(os.path.join(results_dir, "SKIPPED_no_regions.txt")).write_text(
            "No region_annotation column.\n"
        )
        sys.exit(0)

    # Filter out Unlabeled/Bubble regions and Unannotated cell types
    valid_regions = [r for r in adata.obs[region_col].unique()
                     if r not in ("Unlabeled", "Bubble")]
    if len(valid_regions) < 2:
        log.warning("Fewer than 2 valid regions. Skipping.")
        Path(os.path.join(results_dir, "SKIPPED_fewer_than_2_regions.txt")).write_text(
            f"Only {len(valid_regions)} region(s).\n"
        )
        sys.exit(0)

    adata = adata[adata.obs[region_col].isin(valid_regions)].copy()
    log.info("After region filter: %d cells, %d regions: %s",
             adata.n_obs, len(valid_regions), valid_regions)

    # ── Dispatch by analysis level ───────────────────────────────────────
    if analysis_level == "by_region":
        log.info("Level: by_region (all cell types pooled)")
        groupby = [sample_col, region_col]

        counts_df, meta_df = aggregate_pseudobulk(adata, groupby)
        log.info("  Aggregated: %d pseudobulk samples, %d genes",
                 len(counts_df), counts_df.shape[1])

        counts_df, meta_df = filter_pseudobulk(
            counts_df, meta_df, MIN_CELLS, MIN_REPS, region_col
        )
        log.info("  After filtering: %d pseudobulk samples", len(counts_df))

        if len(counts_df) < 4 or meta_df[region_col].nunique() < 2:
            log.warning("  Insufficient data for DE. Skipping.")
            Path(os.path.join(results_dir, "SKIPPED_insufficient.txt")).write_text(
                "Insufficient pseudobulk samples.\n"
            )
            sys.exit(0)

        # Holistic LRT
        run_holistic_lrt(counts_df, meta_df, region_col, results_dir)

        # Pairwise Wald
        run_pairwise_wald(counts_df, meta_df, region_col, results_dir,
                          DE_N_GENES)

    elif analysis_level == "by_celltype_region":
        log.info("Level: by_celltype_region")

        if cell_type_col is None:
            log.warning("No cell type column for this annot_type. Skipping.")
            sys.exit(0)

        # Filter unannotated
        adata = adata[adata.obs[cell_type_col] != "Unannotated"].copy()
        cell_types = [ct for ct in adata.obs[cell_type_col].unique()
                      if ct != "Unannotated"]
        log.info("  Cell types to analyse: %s", cell_types)

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

            groupby = [sample_col, region_col]
            counts_df, meta_df = aggregate_pseudobulk(ct_adata, groupby)
            counts_df, meta_df = filter_pseudobulk(
                counts_df, meta_df, MIN_CELLS, MIN_REPS, region_col
            )

            if len(counts_df) < 4 or meta_df[region_col].nunique() < 2:
                log.warning("    Insufficient pseudobulk samples. Skipping.")
                Path(os.path.join(ct_dir, "SKIPPED_insufficient.txt")).write_text(
                    "Insufficient pseudobulk samples.\n"
                )
                continue

            log.info("    %d pseudobulk samples, %d regions",
                     len(counts_df), meta_df[region_col].nunique())

            run_holistic_lrt(counts_df, meta_df, region_col, ct_dir)
            run_pairwise_wald(counts_df, meta_df, region_col, ct_dir,
                              DE_N_GENES)

    elif analysis_level == "by_niche_region":
        # Future: same logic as by_celltype_region but using niche labels
        log.warning("by_niche_region not yet implemented. "
                     "Run niche identification first.")
        Path(os.path.join(results_dir, "SKIPPED_niche_not_implemented.txt")).write_text(
            "Niche identification must run before this analysis.\n"
        )

    else:
        raise ValueError(f"Unknown analysis_level: {analysis_level}")

    log.info("Pseudobulk DE complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
