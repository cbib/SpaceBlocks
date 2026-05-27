"""
annotate_cells.py – Cell type annotation + scaled-expression refinement
========================================================================
1. TSV-based annotation: maps Leiden clusters → cell types.

2. Scaled-expression refinement for rare/missed cell types:
   - Z-score scales expression ONLY for the marker genes (memory-safe)
   - For each refinement cell type, computes the mean scaled expression
     of its markers per cell
   - A cell is reassigned if:
     (a) mean scaled score > threshold (default 1.0 = 1 SD above mean)
     (b) at least min_markers_expressed genes are individually > 0
         in scaled space
   - Multiple passing → '_mixed'; none passing → keep TSV label
   - Post-filter: types with < min_cells_per_type revert to TSV
   - Scores stored in obs for inspection

3. Generates annotation plots for TSV and refined annotations.
   If cell_type_ingest is present (from ingest_ref), plots it too.
"""

import csv
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
log = logging.getLogger("annotate_cells")


# ── Helpers ──────────────────────────────────────────────────────────────────

def read_tsv_to_dict(tsv_path):
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = {field: [] for field in reader.fieldnames}
        for row in reader:
            for key, val in row.items():
                if val.strip():
                    columns[key].append(val.strip())
    return columns


def find_sample_column(annot_df, sample_id):
    for col in annot_df.columns:
        if col.startswith(sample_id + "_"):
            return col, col[len(sample_id) + 1:]
    if sample_id in annot_df.columns:
        return sample_id, None
    return None, None


def extract_and_scale(adata, gene_names):
    """
    Extract a subset of genes and Z-score scale per gene.
    Returns a dense (n_cells, n_genes) array.  Memory-safe because
    it only densifies the marker columns, not the full matrix.
    """
    idx = [adata.var_names.get_loc(g) for g in gene_names]
    sub = adata.X[:, idx]
    if sparse.issparse(sub):
        sub = sub.toarray()
    else:
        sub = np.asarray(sub, dtype=np.float64)
    means = sub.mean(axis=0)
    stds = sub.std(axis=0)
    stds[stds == 0] = 1
    return (sub - means) / stds


def refine_annotations_scaled(adata, marker_dict, threshold, min_markers,
                              min_cells, tsv_col):
    """
    Refine cell-type annotations using mean scaled expression of markers.
    Only scales the marker genes — not the full expression matrix.
    """
    # Collect all unique marker genes
    all_marker_genes = list({g for genes in marker_dict.values() for g in genes})
    log.info("  Scaling %d unique marker genes (not the full matrix) …", len(all_marker_genes))
    scaled_markers = extract_and_scale(adata, all_marker_genes)

    # Build a gene-name → column-index map for the scaled sub-matrix
    gene_to_col = {g: i for i, g in enumerate(all_marker_genes)}

    # Compute scores per refinement cell type
    score_df = pd.DataFrame(index=adata.obs_names)
    n_positive_df = pd.DataFrame(index=adata.obs_names)

    for ct, genes in marker_dict.items():
        cols = [gene_to_col[g] for g in genes]
        scaled_sub = scaled_markers[:, cols]
        score_df[ct] = scaled_sub.mean(axis=1)
        n_positive_df[ct] = (scaled_sub > 0).sum(axis=1)
        log.info("  %s: score range [%.2f, %.2f], median %.2f",
                 ct, score_df[ct].min(), score_df[ct].max(), score_df[ct].median())

    # Assign
    log.info("  Assigning (threshold=%.2f, min_markers=%d) …", threshold, min_markers)
    refined = []
    for idx in range(adata.n_obs):
        passing = []
        for ct in marker_dict:
            if (score_df[ct].iloc[idx] > threshold
                    and n_positive_df[ct].iloc[idx] >= min_markers):
                passing.append(ct)
        if len(passing) == 1:
            refined.append(passing[0])
        elif len(passing) > 1:
            refined.append("_".join(sorted(passing)) + "_mixed")
        else:
            refined.append(adata.obs[tsv_col].iloc[idx])

    refined = pd.Series(refined, index=adata.obs_names)

    # Post-filter rare new types
    type_counts = refined.value_counts()
    tsv_types = set(adata.obs[tsv_col].astype(str).unique())
    rare_new = [t for t in type_counts[type_counts < min_cells].index
                if t not in tsv_types]
    if rare_new:
        n_rev = refined.isin(rare_new).sum()
        log.info("  Reverting %d cells in %d rare types (< %d cells): %s",
                 n_rev, len(rare_new), min_cells, rare_new)
        mask = refined.isin(rare_new)
        refined[mask] = adata.obs[tsv_col].astype(str)[mask]

    return refined, score_df


def generate_annotation_plots(adata, annot_key, label, plots_dir, sample_id,
                              de_n_genes, library_id=None):
    subdir = os.path.join(plots_dir, label)
    os.makedirs(subdir, exist_ok=True)
    adata.obs[annot_key] = adata.obs[annot_key].astype("category")

    # UMAP
    log.info("  [%s] UMAP …", label)
    sc.pl.umap(adata, color=[annot_key], size=2, wspace=0.25, frameon=False,
               title=f"Cell types ({label})")
    plt.savefig(os.path.join(subdir, f"UMAP_{label}_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # Spatial overview
    try:
        sc.pl.spatial(adata, color=annot_key, spot_size=20, frameon=False,
                      title=f"Spatial – {label}", library_id=library_id)
        plt.savefig(os.path.join(subdir, f"spatial_all_{label}_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        log.warning("  [%s] Spatial failed: %s", label, e)

    # Per-cell-type spatial
    ct_dir = os.path.join(subdir, "spatial_per_celltype")
    os.makedirs(ct_dir, exist_ok=True)
    for ct in adata.obs[annot_key].cat.categories:
        if ct == "Unannotated":
            continue
        adata.obs["_hl"] = adata.obs[annot_key].apply(
            lambda x, c=ct: c if x == c else "Other"
        )
        try:
            sc.pl.spatial(adata, color="_hl", spot_size=20, frameon=False,
                          palette={"Other": "#d3d3d3", ct: "#e41a1c"},
                          title=ct, library_id=library_id)
            safe = ct.replace("/", "_").replace(" ", "_")
            plt.savefig(os.path.join(ct_dir, f"{safe}_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            log.warning("  [%s] Spatial '%s' failed: %s", label, ct, e)
    if "_hl" in adata.obs.columns:
        del adata.obs["_hl"]

    # DE markers
    valid = adata.obs[annot_key].value_counts()
    valid = valid[(valid.index != "Unannotated") & (valid >= 10)]
    if len(valid) >= 2:
        log.info("  [%s] DE (top %d) …", label, de_n_genes)
        sub = adata[adata.obs[annot_key].isin(valid.index)].copy()
        sub.obs[annot_key] = sub.obs[annot_key].cat.remove_unused_categories()
        try:
            sc.tl.rank_genes_groups(sub, groupby=annot_key, method="wilcoxon")
            sc.tl.dendrogram(sub, groupby=annot_key, use_rep="X_pca")

            sc.pl.rank_genes_groups_dotplot(
                sub, groupby=annot_key, standard_scale="var",
                n_genes=de_n_genes, swap_axes=True, dendrogram=True)
            plt.savefig(os.path.join(subdir, f"top{de_n_genes}_dotplot_{label}_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()

            sc.pl.rank_genes_groups_heatmap(
                sub, n_genes=de_n_genes, groupby=annot_key, use_raw=False,
                swap_axes=True, dendrogram=True, show_gene_labels=True)
            plt.savefig(os.path.join(subdir, f"top{de_n_genes}_heatmap_{label}_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()

            sc.pl.rank_genes_groups_matrixplot(
                sub, n_genes=de_n_genes, groupby=annot_key, use_raw=False,
                swap_axes=True, dendrogram=True)
            plt.savefig(os.path.join(subdir, f"top{de_n_genes}_matrixplot_{label}_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()

            df = sc.get.rank_genes_groups_df(sub, None)
            df.to_csv(os.path.join(subdir, f"DE_markers_{label}_{sample_id}.tsv"),
                      sep="\t", index=False)
        except Exception as e:
            log.warning("  [%s] DE failed: %s", label, e)


# ── Parameters ───────────────────────────────────────────────────────────────
sample_id           = snakemake.params.sample_id
threshold           = float(snakemake.params.marker_threshold)
MIN_MARKERS         = int(snakemake.params.min_markers_expressed)
MIN_CELLS_PER_TYPE  = int(snakemake.params.min_cells_per_type)
DE_N_GENES          = int(snakemake.params.de_n_genes)
USE_PRECOMPUTED     = bool(snakemake.params.use_precomputed)
EXT_ANNOT_CFG       = snakemake.params.external_annotation
PRECOMPUTED_DIR     = str(snakemake.params.precomputed_metadata_dir)
ANNOTATION_COLORS   = snakemake.params.annotation_colors
REGION_COLORS       = snakemake.params.region_colors

adata_path     = str(snakemake.input.adata)
metadata_path  = str(snakemake.input.metadata)
annot_tsv_path = str(snakemake.input.cluster_annotations)
markers_path   = str(snakemake.input.annotation_markers)
out_adata_path = str(snakemake.output.adata_annot)
plots_dir      = str(snakemake.output.plots_dir)

try:
    log.info("=" * 70)
    log.info("Annotating sample: %s", sample_id)
    log.info("  Scaled-score threshold: %.2f", threshold)
    log.info("  Min markers expressed:  %d", MIN_MARKERS)
    log.info("  Min cells per type:     %d", MIN_CELLS_PER_TYPE)
    log.info("=" * 70)

    os.makedirs(plots_dir, exist_ok=True)

    adata = sc.read_h5ad(adata_path)
    log.info("Loaded adata: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # Check if ingest labels are present
    has_ingest = "cell_type_ingest" in adata.obs.columns
    if has_ingest:
        log.info("Ingest labels found — will include in plots.")

    library_id = None
    if "spatial" in adata.uns and len(adata.uns["spatial"]) == 1:
        library_id = list(adata.uns["spatial"].keys())[0]

    # Reload precomputed clusters from metadata TSV if configured
    if USE_PRECOMPUTED:
        _ext_meta = os.path.join(PRECOMPUTED_DIR, f"metadata_{sample_id}.tsv") if PRECOMPUTED_DIR else ""
        if _ext_meta and os.path.isfile(_ext_meta):
            _meta_source = _ext_meta
        elif os.path.isfile(metadata_path):
            _meta_source = metadata_path
        else:
            _meta_source = None

        if _meta_source:
            log.info("Reloading clusters from metadata: %s", _meta_source)
            saved = pd.read_csv(_meta_source, sep="\t", index_col=0, comment="#")
            leiden_cols = [c for c in saved.columns if c.startswith("leiden_")]
            for col in leiden_cols:
                if col not in adata.obs.columns:
                    adata.obs[col] = pd.Categorical(saved[col].reindex(adata.obs_names).astype(int))
                    log.info("  Loaded %s from metadata", col)

    # ── 1. TSV-based annotation ──────────────────────────────────────────
    annot_df = pd.read_csv(annot_tsv_path, sep="\t", index_col=0)
    log.info("Annotation TSV columns: %s", list(annot_df.columns))

    annot_col, resolution = find_sample_column(annot_df, sample_id)

    if annot_col is None:
        log.warning("Sample '%s' not found in TSV. All cells → 'Unannotated'.", sample_id)
        adata.obs["cell_type_tsv"] = "Unannotated"
        adata.obs["cell_type_refined"] = "Unannotated"
        Path(out_adata_path).parent.mkdir(parents=True, exist_ok=True)
        adata.write(out_adata_path)
        sys.exit(0)

    log.info("Using column '%s' (resolution: %s)", annot_col,
             resolution if resolution else "not specified")

    if resolution is not None:
        leiden_key = f"leiden_{resolution.replace('.', '_')}"
    else:
        leiden_cols = [c for c in adata.obs.columns if c.startswith("leiden_")]
        if leiden_cols:
            leiden_key = leiden_cols[0]
        else:
            raise ValueError("No leiden columns found in adata.")

    if leiden_key not in adata.obs.columns:
        available = [c for c in adata.obs.columns if c.startswith("leiden_")]
        raise ValueError(
            f"Leiden column '{leiden_key}' not found. Available: {available}."
        )

    cluster_to_celltype = annot_df[annot_col].dropna().to_dict()
    cluster_to_celltype = {str(k): str(v) for k, v in cluster_to_celltype.items()}
    log.info("Cluster → cell type mapping (%d entries): %s",
             len(cluster_to_celltype), cluster_to_celltype)

    adata.obs["cell_type_tsv"] = (
        adata.obs[leiden_key].astype(str)
        .map(cluster_to_celltype)
        .fillna("Unannotated")
        .astype("category")
    )

    tsv_counts = adata.obs["cell_type_tsv"].value_counts()
    log.info("TSV annotation distribution:\n%s", tsv_counts.to_string())

    # ── 2. Marker-based refinement (scaled, memory-safe) ─────────────────
    log.info("Refining with scaled expression scores …")
    raw_marker_dict = read_tsv_to_dict(markers_path)

    for ct, genes in raw_marker_dict.items():
        log.info("  %s: %d marker genes provided", ct, len(genes))

    marker_dict = {}
    for ct, genes in raw_marker_dict.items():
        present = [g for g in genes if g in adata.var_names]
        missing = [g for g in genes if g not in adata.var_names]
        if missing:
            log.warning("  %s: %d/%d NOT FOUND: %s", ct, len(missing), len(genes), missing)
        if present:
            marker_dict[ct] = present
            log.info("  %s: using %d/%d genes", ct, len(present), len(genes))
        else:
            log.warning("  %s: ALL genes missing — excluded", ct)

    if not marker_dict:
        log.warning("No usable markers. Refined = TSV.")
        adata.obs["cell_type_refined"] = adata.obs["cell_type_tsv"].copy()
    else:
        refined_labels, score_df = refine_annotations_scaled(
            adata, marker_dict, threshold, MIN_MARKERS, MIN_CELLS_PER_TYPE,
            tsv_col="cell_type_tsv",
        )
        adata.obs["cell_type_refined"] = pd.Categorical(refined_labels)
        for ct in score_df.columns:
            adata.obs[f"score_{ct}"] = score_df[ct].values

    refined_counts = adata.obs["cell_type_refined"].value_counts()
    log.info("Refined annotation distribution:\n%s", refined_counts.to_string())

    changed = (adata.obs["cell_type_tsv"].astype(str) != adata.obs["cell_type_refined"].astype(str)).sum()
    log.info("  %d cells (%.1f%%) changed by refinement", changed, 100 * changed / adata.n_obs)

    # ── 3. External annotation (from metadata TSV column) ────────────────
    ext_enabled = False
    if isinstance(EXT_ANNOT_CFG, dict) and EXT_ANNOT_CFG.get("enabled", False):
        ext_col = EXT_ANNOT_CFG.get("column", "")
        if ext_col:
            # Resolve metadata source for external annotation
            _ext_meta = os.path.join(PRECOMPUTED_DIR, f"metadata_{sample_id}.tsv") if PRECOMPUTED_DIR else ""
            if _ext_meta and os.path.isfile(_ext_meta):
                _ext_source = _ext_meta
            elif os.path.isfile(metadata_path):
                _ext_source = metadata_path
            else:
                _ext_source = None

            if _ext_source:
                log.info("Loading external annotation from '%s' in %s …", ext_col, _ext_source)
                saved = pd.read_csv(_ext_source, sep="\t", index_col=0, comment="#")
                if ext_col in saved.columns:
                    adata.obs["cell_type_external"] = (
                        saved[ext_col].reindex(adata.obs_names).fillna("Unannotated").astype("category"))
                    ext_enabled = True
                    log.info("  External annotation: %d types",
                             adata.obs["cell_type_external"].nunique())
                else:
                    log.warning("  External column '%s' not in metadata. Skipping.", ext_col)

    # ── 4. Plots ─────────────────────────────────────────────────────────
    log.info("Generating annotation plots …")

    # Apply custom palettes if configured
    if isinstance(ANNOTATION_COLORS, dict):
        for obs_key in ["cell_type_tsv", "cell_type_refined",
                        "cell_type_ingest", "cell_type_external"]:
            cd = ANNOTATION_COLORS.get(obs_key, {})
            if cd and obs_key in adata.obs.columns:
                cats = adata.obs[obs_key].cat.categories
                adata.uns[f"{obs_key}_colors"] = [cd.get(str(c), "#cccccc") for c in cats]

    if REGION_COLORS and "region_annotation" in adata.obs.columns:
        adata.obs["region_annotation"] = adata.obs["region_annotation"].astype("category")
        cats = adata.obs["region_annotation"].cat.categories
        adata.uns["region_annotation_colors"] = [
            REGION_COLORS.get(str(c), "#cccccc") for c in cats
        ]

    generate_annotation_plots(adata, "cell_type_tsv", "tsv", plots_dir,
                              sample_id, DE_N_GENES, library_id)
    generate_annotation_plots(adata, "cell_type_refined", "refined", plots_dir,
                              sample_id, DE_N_GENES, library_id)

    # If ingest labels exist, plot those too
    if has_ingest:
        generate_annotation_plots(adata, "cell_type_ingest", "ingest", plots_dir,
                                  sample_id, DE_N_GENES, library_id)

    if ext_enabled:
        generate_annotation_plots(adata, "cell_type_external", "external", plots_dir,
                                  sample_id, DE_N_GENES, library_id)

    # Side-by-side comparison (dynamic number of panels)
    annot_cols_present = [c for c in ["cell_type_tsv", "cell_type_refined",
                                       "cell_type_ingest", "cell_type_external"]
                          if c in adata.obs.columns]
    n_panels = len(annot_cols_present)
    if n_panels >= 2:
        fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 6))
        for ax, col in zip(axes, annot_cols_present):
            sc.pl.umap(adata, color=col, size=2, frameon=False,
                       title=col.replace("cell_type_", ""), ax=ax, show=False)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, f"UMAP_comparison_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    # Side-by-side SPATIAL comparison (only if >1 annotation method)
    if n_panels >= 2:
        try:
            fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 6))
            for ax, col in zip(axes, annot_cols_present):
                sc.pl.spatial(adata, color=col, spot_size=20, frameon=False,
                              title=col.replace("cell_type_", ""),
                              library_id=library_id, ax=ax, show=False)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"spatial_comparison_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            log.warning("Spatial comparison plot failed: %s", e)

    # Per-sample composition barplots (by region)
    has_regions = ("region_annotation" in adata.obs.columns
                   and adata.obs["region_annotation"].nunique() > 1
                   and not all(adata.obs["region_annotation"] == "Unlabeled"))
    if has_regions:
        for annot_col in annot_cols_present:
            try:
                ct = pd.crosstab(adata.obs["region_annotation"], adata.obs[annot_col])
                ct_norm = ct.div(ct.sum(axis=1), axis=0) * 100
                label = annot_col.replace("cell_type_", "")

                # Apply custom colors if available
                colors = None
                if isinstance(ANNOTATION_COLORS, dict):
                    cd = ANNOTATION_COLORS.get(annot_col, {})
                    if cd:
                        colors = [cd.get(str(c), "#cccccc") for c in ct.columns]

                for data, suffix, ylabel in [(ct, "absolute", "Number of cells"),
                                              (ct_norm, "relative", "Percentage (%)")]:
                    ax = data.plot(kind="bar", stacked=True,
                                   figsize=(max(8, len(ct) * 1.2), 6),
                                   color=colors)
                    ax.set_ylabel(ylabel)
                    ax.set_xlabel("Region")
                    ax.legend(title=label, bbox_to_anchor=(1.05, 1), loc="upper left")
                    plt.tight_layout()
                    plt.savefig(os.path.join(plots_dir,
                                f"barplot_{label}_{suffix}_{sample_id}.png"),
                                dpi=300, bbox_inches="tight")
                    plt.close()
            except Exception as e:
                log.warning("Barplot for %s failed: %s", annot_col, e)

    # Score distributions
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if score_cols:
        fig, axes = plt.subplots(1, len(score_cols),
                                 figsize=(5 * len(score_cols), 4))
        if len(score_cols) == 1:
            axes = [axes]
        for ax, col in zip(axes, score_cols):
            ax.hist(adata.obs[col].values, bins=50, edgecolor="black", alpha=0.7)
            ax.axvline(threshold, color="red", ls="--", label=f"threshold={threshold}")
            ax.set_title(col.replace("score_", ""))
            ax.set_xlabel("Mean scaled score")
            ax.set_ylabel("Cells")
            ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, f"score_distributions_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    # ── 5. Save ──────────────────────────────────────────────────────────
    log.info("Saving annotated adata → %s", out_adata_path)
    Path(out_adata_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_adata_path)

    log.info("Annotation complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
