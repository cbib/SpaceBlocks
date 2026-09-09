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
import re
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
    "external_annotation": "cell_type_external",
    "refined_annotation": "cell_type_refined",
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


def _safe_name(value):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return safe or "group"


def _group_label(row, columns):
    if len(columns) == 1:
        return str(row[columns[0]])
    return " | ".join(f"{column}={row[column]}" for column in columns)


class DesignMetadataError(ValueError):
    """Requested model metadata is not constant within a count unit."""


def _attach_design_metadata(pdata, source, sample_col, groups_col, columns):
    """Attach metadata that is constant within every pseudobulk count unit.

    decoupler retains grouping columns, but preservation of unrelated obs columns
    is version-dependent. Reconstructing them explicitly makes the DE contract
    stable and detects invalid designs where a requested value varies inside one
    pseudobulk observation.
    """
    columns = list(dict.fromkeys(c for c in columns if c not in {sample_col, groups_col}))
    if not columns:
        return pdata

    keys = [sample_col, groups_col]
    frame = source.obs[keys + columns].copy()
    for column in keys + columns:
        frame[column] = frame[column].astype(str)

    grouped = frame.groupby(keys, observed=True, dropna=False)
    for column in columns:
        variable = grouped[column].nunique(dropna=False)
        if (variable > 1).any():
            examples = [" / ".join(map(str, idx if isinstance(idx, tuple) else [idx]))
                        for idx in variable[variable > 1].index[:5]]
            raise DesignMetadataError(
                f"Design column '{column}' varies within pseudobulk unit(s): {examples}. "
                "Use an aggregation that keeps this variable separate."
            )

    lookup = grouped[columns].first()
    pdata_keys = pdata.obs[keys].astype(str)
    lookup_index = pd.MultiIndex.from_frame(pdata_keys)
    for column in columns:
        pdata.obs[column] = lookup[column].reindex(lookup_index).to_numpy()
    return pdata


def _aggregate(source, groups_col, design_columns):
    pdata = dc.get_pseudobulk(
        source,
        sample_col="sample",
        groups_col=groups_col,
        layer="raw_counts",
        mode="sum",
        min_cells=MIN_CELLS,
        min_counts=MIN_COUNTS,
    )
    return _attach_design_metadata(
        pdata, source, "sample", groups_col, design_columns
    )


# ── Parameters ───────────────────────────────────────────────────────────────
annot_type = str(snakemake.params.annot_type)
analysis_spec = dict(snakemake.params.analysis_spec)
analysis_name = str(analysis_spec["name"])
aggregation = str(analysis_spec["aggregation"])
group_by = list(analysis_spec["group_by"])
paired_by = str(analysis_spec.get("paired_by", "") or "")
covariates = list(analysis_spec.get("covariates") or [])
exclude_levels = dict(analysis_spec.get("exclude_levels") or {})
MIN_CELLS = int(snakemake.params.min_cells_per_pseudobulk)
MIN_COUNTS = int(snakemake.params.min_counts_per_pseudobulk)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS = snakemake.params.region_colors
DPI = int(getattr(snakemake.params, "dpi", 300))
EXTRA_ANNOT_COLUMNS = list(getattr(snakemake.params, "extra_annot_columns", []) or [])
SAMPLE_COLORS = getattr(snakemake.params, "sample_colors", {}) or {}
agg_dir = str(snakemake.output.agg_dir)

try:
    log.info("=" * 70)
    log.info("Pseudobulk aggregation: %s / %s", annot_type, analysis_name)
    log.info("  aggregation=%s, group_by=%s", aggregation, group_by)
    log.info("  min_cells=%d, min_counts=%d", MIN_CELLS, MIN_COUNTS)
    log.info("=" * 70)

    matrices_dir = os.path.join(agg_dir, "matrices")
    metadata_dir = os.path.join(agg_dir, "metadata")
    plots_dir = os.path.join(agg_dir, "plots")
    for directory in [agg_dir, matrices_dir, metadata_dir, plots_dir]:
        os.makedirs(directory, exist_ok=True)

    adata = sc.read_h5ad(str(snakemake.input.adata))
    log.info("Loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # sample_batch is created by integrate_samples from the authoritative core
    # sample sheet. Prefer it over a potentially renamed contract obs['sample'].
    if "sample_batch" in adata.obs.columns:
        adata.obs["sample"] = adata.obs["sample_batch"].astype(str).to_numpy()
    elif "sample" not in adata.obs.columns:
        raise ValueError("No canonical sample column found.")

    cell_type_col = ANNOT_COL_MAP.get(annot_type)
    referenced_columns = set(group_by + covariates + list(exclude_levels))
    if paired_by:
        referenced_columns.add(paired_by)
    needs_cell_type = (
        aggregation in {"by_celltype", "by_celltype_region"}
        or "cell_type" in referenced_columns
    )
    if needs_cell_type:
        if not cell_type_col or cell_type_col not in adata.obs.columns:
            log.warning("Column '%s' not found. Skipping.", cell_type_col)
            Path(os.path.join(agg_dir, f"SKIPPED_no_{cell_type_col}.txt")).write_text(
                f"Column {cell_type_col} not found.\n"
            )
            sys.exit(0)
        adata.obs["cell_type"] = adata.obs[cell_type_col].astype(str).to_numpy()

    needs_region = aggregation in {"by_region", "by_celltype_region"}
    if needs_region and "region_annotation" not in adata.obs.columns:
        log.warning("No region_annotation for region-based aggregation. Skipping.")
        Path(os.path.join(agg_dir, "SKIPPED_no_regions.txt")).write_text("")
        sys.exit(0)

    required_columns = list(dict.fromkeys(
        ["sample"] + group_by + covariates + ([paired_by] if paired_by else [])
    ))
    missing = [column for column in required_columns if column not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Pseudobulk analysis '{analysis_name}' requires missing obs column(s): {missing}. "
            "Add them to core_samples.tsv or correct group_by/covariates/paired_by."
        )

    for column in required_columns:
        adata = _clean_nan(adata, column)

    # Exclusions are shared by aggregation, explicit Wald contrasts, and LRT.
    for column, excluded in exclude_levels.items():
        if column not in adata.obs.columns:
            raise ValueError(
                f"exclude_levels references missing column '{column}' in analysis '{analysis_name}'."
            )
        excluded = {str(value) for value in (excluded or [])}
        if excluded:
            values = adata.obs[column].astype(str)
            keep = ~values.isin(excluded)
            log.info("  Excluding %d cells from %s levels: %s",
                     int((~keep).sum()), column, sorted(excluded))
            adata = adata[keep].copy()

    if adata.n_obs == 0:
        Path(os.path.join(agg_dir, "SKIPPED_no_cells_after_exclusion.txt")).write_text("")
        sys.exit(0)
    if "raw_counts" not in adata.layers:
        raise ValueError("raw_counts layer missing.")

    group_palette = {}
    if len(group_by) == 1:
        group_column = group_by[0]
        if group_column == "region_annotation":
            group_palette = REGION_COLORS if isinstance(REGION_COLORS, dict) else {}
        elif isinstance(SAMPLE_COLORS, dict):
            group_palette = SAMPLE_COLORS.get(group_column, {}) or {}

    # Plot-only annotations are optional. Model columns were validated above and
    # remain mandatory, while absent plotting annotations should not abort DE.
    present_extra_annotations = [
        column for column in EXTRA_ANNOT_COLUMNS if column in adata.obs.columns
    ]
    design_columns = list(dict.fromkeys(
        group_by + covariates + ([paired_by] if paired_by else []) + present_extra_annotations
    ))
    manifest_rows = []
    used_prefixes = set()

    if aggregation in {"all_cells", "by_celltype"}:
        adata.obs["__all__"] = "all"
        groups_col = "__all__"
    else:
        groups_col = "region_annotation"

    if aggregation in {"all_cells", "by_region"}:
        subgroups = [("pooled", "all_cells", adata)]
    else:
        subgroups = []
        for cell_type in sorted(adata.obs["cell_type"].astype(str).unique()):
            prefix = _safe_name(cell_type)
            if prefix in used_prefixes:
                raise ValueError(
                    f"Cell-type names collide after filename sanitization: '{cell_type}' -> '{prefix}'."
                )
            used_prefixes.add(prefix)
            subset = adata[adata.obs["cell_type"].astype(str) == cell_type].copy()
            subgroups.append((prefix, cell_type, subset))

    for prefix, grouping, subset in subgroups:
        if subset.n_obs < MIN_CELLS:
            log.warning("  '%s': too few cells (%d). Skipping.", grouping, subset.n_obs)
            continue
        try:
            pdata = _aggregate(subset, groups_col, design_columns)
        except DesignMetadataError:
            raise
        except Exception as exc:
            log.warning("  '%s': aggregation failed: %s", grouping, exc)
            continue
        if pdata.n_obs == 0:
            log.warning("  '%s': no pseudobulk observations passed QC. Skipping.", grouping)
            continue

        pdata.obs[".comparison_group"] = [
            _group_label(row, group_by) for _, row in pdata.obs.iterrows()
        ]
        _save_pseudobulk(pdata, matrices_dir, metadata_dir, prefix)
        _plot_pseudobulk_qc(
            pdata, ".comparison_group", "sample", prefix, plots_dir,
            group_palette, ANNOTATION_COLORS,
            EXTRA_ANNOT_COLUMNS, SAMPLE_COLORS,
        )
        manifest_rows.append({
            "prefix": prefix,
            "grouping": grouping,
            "n_samples": pdata.n_obs,
            "n_genes": pdata.n_vars,
            "condition_col": ".comparison_group",
            "sample_col": "sample",
            "analysis_name": analysis_name,
            "aggregation": aggregation,
        })

    if manifest_rows:
        pd.DataFrame(manifest_rows).to_csv(
            os.path.join(agg_dir, "manifest.tsv"), sep="\t", index=False
        )
    else:
        Path(os.path.join(agg_dir, "SKIPPED_no_pseudobulks.txt")).write_text("")

    log.info("Aggregation complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
