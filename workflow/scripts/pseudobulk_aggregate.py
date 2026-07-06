"""
pseudobulk_aggregate.py – Aggregate cells into pseudobulk matrices
====================================================================
Output structure:
  aggregated/
  ├── manifest.tsv
  ├── matrices/     count TSVs
  ├── metadata/     metadata TSVs
  └── plots/        3-panel QC per subgroup
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
from scipy import sparse
from sklearn.decomposition import PCA
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


def _save_pseudobulk(pdata, matrices_dir, metadata_dir, prefix):
    if sparse.issparse(pdata.X):
        counts = pd.DataFrame(pdata.X.toarray(), index=pdata.obs_names,
                              columns=pdata.var_names)
    else:
        counts = pd.DataFrame(np.asarray(pdata.X), index=pdata.obs_names,
                              columns=pdata.var_names)
    counts = counts.round().astype(int)
    counts.to_csv(os.path.join(matrices_dir, f"counts_{prefix}.tsv"), sep="\t")
    pdata.obs.to_csv(os.path.join(metadata_dir, f"metadata_{prefix}.tsv"), sep="\t")
    log.info("  Saved: %s (%d samples × %d genes)", prefix, pdata.n_obs, pdata.n_vars)


def _plot_pseudobulk_qc(pdata, condition_col, sample_col, prefix, plots_dir,
                        region_colors=None, annotation_colors=None,
                        extra_annot_columns=None, sample_colors=None):
    """
    QC panels: dc.plot_psbulk_samples + PCA by condition + PCA by sample, plus one
    extra PCA panel per design column (coloured from sample_colors, grey fallback).
    Legends placed outside the plot area. Config palettes (region_colors for the
    condition PCA, annotation_colors['sample_batch'] for the sample PCA) are used
    when available, with colormap fallbacks otherwise.
    """
    design_present = [c for c in (extra_annot_columns or []) if c in pdata.obs.columns]
    sample_colors = sample_colors if isinstance(sample_colors, dict) else {}
    n_panels = 3 + len(design_present)
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels + 1, 5))
    fig.suptitle(f"Pseudobulk QC — {prefix}", fontsize=14, fontweight="bold")

    # Panel 1: cells-vs-counts per pseudobulk sample, coloured by condition
    try:
        rc = region_colors if isinstance(region_colors, dict) else {}
        _conds = sorted(set(pdata.obs[condition_col].astype(str)))
        _cmap_c = plt.cm.get_cmap("Set1", max(len(_conds), 1))
        p1_colors = {c: rc.get(c, _cmap_c(i)) for i, c in enumerate(_conds)}

        if {"psbulk_n_cells", "psbulk_counts"}.issubset(pdata.obs.columns):
            cvals = pdata.obs[condition_col].astype(str)
            for c in _conds:
                m = (cvals == c).values
                axes[0].scatter(pdata.obs.loc[m, "psbulk_n_cells"],
                                pdata.obs.loc[m, "psbulk_counts"],
                                c=[p1_colors[c]], label=c, s=50,
                                edgecolors="black", linewidths=0.5)
            axes[0].set_xscale("log")
            axes[0].set_yscale("log")
            axes[0].set_xlabel("Number of cells")
            axes[0].set_ylabel("Total counts")
            axes[0].legend(fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5),
                           borderaxespad=0, frameon=False)
        else:
            # QC columns absent — fall back to decoupler's own plot
            dc.plot_psbulk_samples(pdata, groupby=condition_col, ax=axes[0])
            axes[0].tick_params(axis="x", rotation=45)
        axes[0].set_title("Cells per pseudobulk sample")
    except Exception as e:
        log.warning("  plot_psbulk_samples failed: %s", e)
        axes[0].text(0.5, 0.5, f"Plot failed:\n{e}", ha="center", va="center",
                     transform=axes[0].transAxes, fontsize=8)

    # PCA computation
    try:
        pdata_norm = pdata.copy()
        sc.pp.normalize_total(pdata_norm, target_sum=1e6)
        sc.pp.log1p(pdata_norm)

        if sparse.issparse(pdata_norm.X):
            X_dense = pdata_norm.X.toarray()
        else:
            X_dense = np.asarray(pdata_norm.X)

        gene_var = X_dense.var(axis=0)
        X_dense = X_dense[:, gene_var > 0]

        n_components = min(2, X_dense.shape[0], X_dense.shape[1])
        if n_components < 2:
            raise ValueError(f"Too few samples/genes for PCA ({X_dense.shape})")

        pca = PCA(n_components=n_components)
        pcs = pca.fit_transform(X_dense)
        var_explained = pca.explained_variance_ratio_ * 100

        conditions = pdata.obs[condition_col].values
        samples = pdata.obs[sample_col].values if sample_col in pdata.obs.columns else None

        # Panel 2: PCA by condition
        unique_conds = sorted(set(conditions))
        rc = region_colors if isinstance(region_colors, dict) else {}
        _cmap = plt.cm.get_cmap("Set1", max(len(unique_conds), 1))
        cond_colors = {c: rc.get(str(c), _cmap(i))
                       for i, c in enumerate(unique_conds)}

        for cond in unique_conds:
            mask = conditions == cond
            axes[1].scatter(pcs[mask, 0], pcs[mask, 1], c=[cond_colors[cond]],
                            label=cond, s=60, edgecolors="black", linewidths=0.5)
        axes[1].set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
        axes[1].set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
        axes[1].set_title(f"PCA — {condition_col}")
        axes[1].legend(fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5),
                       borderaxespad=0, frameon=False)

        # Panel 3: PCA by sample
        if samples is not None:
            unique_samples = sorted(set(samples))
            sc_cd = (annotation_colors.get("sample_batch", {})
                     if isinstance(annotation_colors, dict) else {})
            cmap_s = plt.cm.get_cmap("tab20", max(len(unique_samples), 1))
            sample_pal = {s: sc_cd.get(str(s), cmap_s(i))
                          for i, s in enumerate(unique_samples)}
            for s in unique_samples:
                mask = samples == s
                axes[2].scatter(pcs[mask, 0], pcs[mask, 1], c=[sample_pal[s]],
                                label=s, s=60, edgecolors="black", linewidths=0.5)
            axes[2].set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
            axes[2].set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
            axes[2].set_title(f"PCA — {sample_col}")
            axes[2].legend(fontsize=5, loc="center left", bbox_to_anchor=(1.02, 0.5),
                           borderaxespad=0, frameon=False, ncol=1)
        else:
            axes[2].text(0.5, 0.5, "No sample column", ha="center", va="center",
                         transform=axes[2].transAxes)

        # Extra PCA panels: one per design column (grey for palette-less values)
        for _i, _dc in enumerate(design_present):
            ax = axes[3 + _i]
            dvals = pdata.obs[_dc].astype(str).values
            duniq = sorted(set(dvals))
            dpal = sample_colors.get(_dc, {}) if isinstance(sample_colors, dict) else {}
            dcolors = {v: dpal.get(str(v), "#cccccc") for v in duniq}
            for v in duniq:
                m = dvals == v
                ax.scatter(pcs[m, 0], pcs[m, 1], c=[dcolors[v]], label=v, s=60,
                           edgecolors="black", linewidths=0.5)
            ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
            ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
            ax.set_title(f"PCA — {_dc}")
            ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5),
                      borderaxespad=0, frameon=False)

    except Exception as e:
        log.warning("  PCA plots failed: %s", e)
        for ax_idx in range(1, n_panels):
            axes[ax_idx].text(0.5, 0.5, f"PCA failed:\n{e}", ha="center",
                              va="center", transform=axes[ax_idx].transAxes, fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    safe_prefix = prefix.replace("/", "_").replace(" ", "_")
    plt.savefig(os.path.join(plots_dir, f"{safe_prefix}_qc.png"),
                dpi=DPI, bbox_inches="tight")
    plt.close()


# ── Parameters ───────────────────────────────────────────────────────────────
annot_type     = snakemake.params.annot_type
analysis_level = snakemake.params.analysis_level
MIN_CELLS      = int(snakemake.params.min_cells_per_pseudobulk)
MIN_COUNTS     = int(snakemake.params.min_counts_per_pseudobulk)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors
DPI          = int(getattr(snakemake.params, "dpi", 300))
EXTRA_ANNOT_COLUMNS = list(getattr(snakemake.params, "extra_annot_columns", []) or [])
SAMPLE_COLORS  = getattr(snakemake.params, "sample_colors", {}) or {}
agg_dir        = str(snakemake.output.agg_dir)

try:
    log.info("=" * 70)
    log.info("Pseudobulk aggregation: %s / %s", annot_type, analysis_level)
    log.info("  min_cells=%d, min_counts=%d", MIN_CELLS, MIN_COUNTS)
    log.info("=" * 70)

    matrices_dir = os.path.join(agg_dir, "matrices")
    metadata_dir = os.path.join(agg_dir, "metadata")
    plots_dir    = os.path.join(agg_dir, "plots")
    for d in [agg_dir, matrices_dir, metadata_dir, plots_dir]:
        os.makedirs(d, exist_ok=True)

    adata = sc.read_h5ad(str(snakemake.input.adata))
    log.info("Loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)

    cell_type_col = ANNOT_COL_MAP.get(annot_type)
    if cell_type_col and cell_type_col not in adata.obs.columns:
        log.warning("Column '%s' not found. Skipping.", cell_type_col)
        Path(os.path.join(agg_dir, f"SKIPPED_no_{cell_type_col}.txt")).write_text(
            f"Column {cell_type_col} not found.\n"
        )
        sys.exit(0)

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

    for col in [region_col, sample_col]:
        adata = _clean_nan(adata, col)
    if cell_type_col:
        adata = _clean_nan(adata, cell_type_col)

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
        _save_pseudobulk(pdata, matrices_dir, metadata_dir, "pooled")
        _plot_pseudobulk_qc(pdata, region_col, sample_col, "pooled", plots_dir,
                            REGION_COLORS, ANNOTATION_COLORS,
                            EXTRA_ANNOT_COLUMNS, SAMPLE_COLORS)
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

            if ct_adata.n_obs < MIN_CELLS:
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

            _save_pseudobulk(pdata, matrices_dir, metadata_dir, safe_ct)
            _plot_pseudobulk_qc(pdata, region_col, sample_col, safe_ct, plots_dir,
                                REGION_COLORS, ANNOTATION_COLORS,
                                EXTRA_ANNOT_COLUMNS, SAMPLE_COLORS)
            manifest_rows.append({
                "prefix": safe_ct, "grouping": ct,
                "n_samples": pdata.n_obs, "n_genes": pdata.n_vars,
                "condition_col": region_col, "sample_col": sample_col,
            })

    elif analysis_level == "by_niche_region":
        log.warning("by_niche_region not yet implemented.")
        Path(os.path.join(agg_dir, "SKIPPED_niche_not_implemented.txt")).write_text("")

    if manifest_rows:
        pd.DataFrame(manifest_rows).to_csv(
            os.path.join(agg_dir, "manifest.tsv"), sep="\t", index=False
        )

    log.info("Aggregation complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
