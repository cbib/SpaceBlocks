"""
annotate_cells.py – Cell type annotation
=========================================
1. TSV-based annotation: maps Leiden clusters → cell types (cell_type_tsv).

2. External annotation (optional): loads a cell-type column from the metadata
   TSV (cell_type_external) when enabled in config.

3. Generates annotation plots for the TSV annotation. If cell_type_ingest
   (from ingest_ref) or cell_type_external are present, plots them too, plus
   side-by-side UMAP/spatial comparisons and per-sample composition barplots.
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

# Shared composition-barplot helpers (black edge, legend==stack order, config colours)
try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:                      # very old Snakemake
    _here = os.getcwd()
sys.path.insert(0, _here)
from composition_barplots import composition_pair, find_niche_column


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
                dpi=DPI, bbox_inches="tight")
    plt.close()

    # Spatial overview
    try:
        sc.pl.spatial(adata, color=annot_key, spot_size=20, frameon=False,
                      title=f"Spatial – {label}", library_id=library_id)
        plt.savefig(os.path.join(subdir, f"spatial_all_{label}_{sample_id}.png"),
                    dpi=DPI, bbox_inches="tight")
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
                          palette={"Other": "#d3d3d3", ct: "#000000"},
                          title=ct, library_id=library_id)
            safe = ct.replace("/", "_").replace(" ", "_")
            plt.savefig(os.path.join(ct_dir, f"{safe}_{sample_id}.png"),
                        dpi=DPI, bbox_inches="tight")
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
                        dpi=DPI, bbox_inches="tight")
            plt.close()

            sc.pl.rank_genes_groups_heatmap(
                sub, n_genes=de_n_genes, groupby=annot_key, use_raw=False,
                swap_axes=True, dendrogram=True, show_gene_labels=True)
            plt.savefig(os.path.join(subdir, f"top{de_n_genes}_heatmap_{label}_{sample_id}.png"),
                        dpi=DPI, bbox_inches="tight")
            plt.close()

            sc.pl.rank_genes_groups_matrixplot(
                sub, n_genes=de_n_genes, groupby=annot_key, use_raw=False,
                swap_axes=True, dendrogram=True)
            plt.savefig(os.path.join(subdir, f"top{de_n_genes}_matrixplot_{label}_{sample_id}.png"),
                        dpi=DPI, bbox_inches="tight")
            plt.close()

            df = sc.get.rank_genes_groups_df(sub, None)
            df.to_csv(os.path.join(subdir, f"DE_markers_{label}_{sample_id}.tsv"),
                      sep="\t", index=False)
        except Exception as e:
            log.warning("  [%s] DE failed: %s", label, e)


# ── Parameters ───────────────────────────────────────────────────────────────
sample_id           = snakemake.params.sample_id
MIN_CELLS_PER_TYPE  = int(snakemake.params.min_cells_per_type)
DE_N_GENES          = int(snakemake.params.de_n_genes)
USE_PRECOMPUTED     = bool(snakemake.params.use_precomputed)
EXT_ANNOT_CFG       = snakemake.params.external_annotation
PRECOMPUTED_DIR     = str(snakemake.params.precomputed_metadata_dir)
ANNOTATION_COLORS   = snakemake.params.annotation_colors
REGION_COLORS       = snakemake.params.region_colors
DPI          = int(getattr(snakemake.params, "dpi", 300))
NICHE_COLUMN        = getattr(snakemake.params, "niche_column", "")

adata_path     = str(snakemake.input.adata)
metadata_path  = str(snakemake.input.metadata)
annot_tsv_path = str(snakemake.input.cluster_annotations)
out_adata_path = str(snakemake.output.adata_annot)
plots_dir      = str(snakemake.output.plots_dir)

try:
    log.info("=" * 70)
    log.info("Annotating sample: %s", sample_id)
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

    # ── 2. External annotation (from metadata TSV column) ────────────────
    # When enabled the pipeline MUST find the column for this sample — a missing
    # source or column is a hard error (fail loudly rather than silently skipping),
    # so a run only proceeds when the external labels are in place for every sample.
    ext_enabled = False
    if isinstance(EXT_ANNOT_CFG, dict) and EXT_ANNOT_CFG.get("enabled", False):
        ext_col = EXT_ANNOT_CFG.get("column", "")
        if not ext_col:
            raise ValueError("external_annotation.enabled is true but no 'column' is set in config.")
        # Resolve metadata source for external annotation
        _ext_meta = os.path.join(PRECOMPUTED_DIR, f"metadata_{sample_id}.tsv") if PRECOMPUTED_DIR else ""
        if _ext_meta and os.path.isfile(_ext_meta):
            _ext_source = _ext_meta
        elif os.path.isfile(metadata_path):
            _ext_source = metadata_path
        else:
            _ext_source = None

        if _ext_source is None:
            raise FileNotFoundError(
                f"external_annotation enabled but no metadata found for sample "
                f"'{sample_id}' (looked in precomputed_metadata_dir and {metadata_path}).")

        log.info("Loading external annotation from '%s' in %s …", ext_col, _ext_source)
        saved = pd.read_csv(_ext_source, sep="\t", index_col=0, comment="#")
        if ext_col not in saved.columns:
            raise ValueError(
                f"external_annotation column '{ext_col}' not found in {_ext_source} "
                f"(sample '{sample_id}'). Available columns: {list(saved.columns)}")
        adata.obs["cell_type_external"] = (
            saved[ext_col].reindex(adata.obs_names).fillna("Unannotated").astype("category"))
        ext_enabled = True
        log.info("  External annotation: %d types (%d cells unmatched → Unannotated)",
                 adata.obs["cell_type_external"].nunique(),
                 int((adata.obs["cell_type_external"] == "Unannotated").sum()))

    # ── 2b. Spatial niche labels (from spatial_niches rule / external) ───
    _ni = getattr(snakemake.input, "spatial_niche", "")
    if isinstance(_ni, (list, tuple)):
        niche_tsv = str(_ni[0]) if len(_ni) else ""
    else:
        niche_tsv = str(_ni) if _ni else ""
    if niche_tsv and os.path.isfile(niche_tsv):
        log.info("Merging spatial niche labels from %s …", niche_tsv)
        ndf = pd.read_csv(niche_tsv, sep="\t", index_col=0)
        col = "spatial_niche" if "spatial_niche" in ndf.columns else ndf.columns[0]
        mapped = ndf[col].reindex(adata.obs_names)
        n_assigned = int(mapped.notna().sum())
        mapped = mapped.fillna("Unassigned").astype(str)
        adata.obs["spatial_niche"] = pd.Categorical(mapped)
        log.info("  spatial_niche: %d/%d cells assigned, %d niches",
                 n_assigned, adata.n_obs, adata.obs["spatial_niche"].nunique())

    # ── 3. Plots ─────────────────────────────────────────────────────────
    log.info("Generating annotation plots …")

    # Apply custom palettes if configured
    if isinstance(ANNOTATION_COLORS, dict):
        for obs_key in ["cell_type_tsv", "cell_type_ingest", "cell_type_external"]:
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

    # If ingest labels exist, plot those too
    if has_ingest:
        generate_annotation_plots(adata, "cell_type_ingest", "ingest", plots_dir,
                                  sample_id, DE_N_GENES, library_id)

    if ext_enabled:
        generate_annotation_plots(adata, "cell_type_external", "external", plots_dir,
                                  sample_id, DE_N_GENES, library_id)

    # Side-by-side comparison (dynamic number of panels)
    annot_cols_present = [c for c in ["cell_type_tsv", "cell_type_ingest",
                                       "cell_type_external"]
                          if c in adata.obs.columns]
    n_panels = len(annot_cols_present)
    if n_panels >= 2:
        fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 6),
                                 gridspec_kw={"wspace": 0.5})
        for ax, col in zip(axes, annot_cols_present):
            sc.pl.umap(adata, color=col, size=2, frameon=False,
                       title=col.replace("cell_type_", ""), ax=ax, show=False,
                       legend_fontsize=6, na_in_legend=False)
        plt.savefig(os.path.join(plots_dir, f"UMAP_comparison_{sample_id}.png"),
                    dpi=DPI, bbox_inches="tight")
        plt.close()

    # Side-by-side SPATIAL comparison (only if >1 annotation method)
    if n_panels >= 2:
        try:
            fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 6),
                                     gridspec_kw={"wspace": 0.5})
            for ax, col in zip(axes, annot_cols_present):
                sc.pl.spatial(adata, color=col, spot_size=20, frameon=False,
                              title=col.replace("cell_type_", ""),
                              library_id=library_id, ax=ax, show=False,
                              legend_fontsize=6, na_in_legend=False)
            plt.savefig(os.path.join(plots_dir, f"spatial_comparison_{sample_id}.png"),
                        dpi=DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            log.warning("Spatial comparison plot failed: %s", e)

    # ── Per-sample composition barplots (sample / region / niche) ────────
    annot_cols_present = [c for c in ["cell_type_tsv", "cell_type_ingest",
                                      "cell_type_external"]
                          if c in adata.obs.columns]

    # sample column (one bar per sample; here a single sample). Per-sample adatas
    # carry "sample"; "sample_batch" is an integration-time concat label and does
    # not exist here, so it is not consulted.
    sample_col = "sample" if "sample" in adata.obs.columns else None
    if sample_col is None:
        adata.obs["sample"] = sample_id
        sample_col = "sample"

    has_regions = ("region_annotation" in adata.obs.columns
                   and adata.obs["region_annotation"].nunique() > 1
                   and not all(adata.obs["region_annotation"] == "Unlabeled"))
    niche_col = find_niche_column(adata, NICHE_COLUMN)   # NICHE_COLUMN may be ""

    bar_dir = os.path.join(plots_dir, "composition")
    for annot_col in annot_cols_present:
        cmap = (ANNOTATION_COLORS.get(annot_col, {})
                if isinstance(ANNOTATION_COLORS, dict) else {})
        label = annot_col.replace("cell_type_", "")
        composition_pair(adata, sample_col, annot_col, cmap, bar_dir,
                         f"{label}_by_sample", group_label="Sample",
                         cat_label=label, dpi=DPI)
        if has_regions:
            composition_pair(adata, "region_annotation", annot_col, cmap,
                             bar_dir, f"{label}_by_region",
                             group_label="Region", cat_label=label, dpi=DPI)
        if niche_col:
            composition_pair(adata, niche_col, annot_col, cmap, bar_dir,
                             f"{label}_by_niche", group_label="Niche",
                             cat_label=label, dpi=DPI)

    # ── 4. Save ──────────────────────────────────────────────────────────
    log.info("Saving annotated adata → %s", out_adata_path)
    adata.uns["annotation_leiden_key"] = leiden_key   # the TSV-annotation res
    Path(out_adata_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_adata_path)

    log.info("Annotation complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
