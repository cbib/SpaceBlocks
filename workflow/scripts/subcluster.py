"""
subcluster.py – Subset and subcluster a cell compartment
==========================================================
Output per branch (Harmony / NoHarmony):
  ├── adata_{subcompartment}_{branch}.h5ad
  ├── cell_counts_summary.tsv
  ├── Clustering_evaluation/   silhouette + clustree
  ├── QC_plots/                QC UMAPs + violins per cluster
  ├── UMAPs/                   cluster, sample, region, split UMAPs
  ├── ClusterMarkers/          top-N DE dotplot, heatmap, matrixplot, TSV
  └── Barplots/                stacked barplots (abs + rel) for all combos
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
from pyclustree import clustree
import scanpy as sc
import scanpy.external as sce
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.preprocessing import StandardScaler

# Shared composition-barplot helper (black edge, legend==stack order, config colours)
try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:                      # very old Snakemake
    _here = os.getcwd()
sys.path.insert(0, _here)
from composition_barplots import save_stacked_composition


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
log = logging.getLogger("subcluster")

RANDOM_SEED = int(snakemake.params.random_seed)

# ── Helpers ──────────────────────────────────────────────────────────────────

def _resolution_range(rmin, rmax, step):
    n = round((rmax - rmin) / step)
    return [round(rmin + i * step, 10) for i in range(n + 1)]


def split_umap(adata, split_by, ncol=None, nrow=None, **kwargs):
    categories = adata.obs[split_by].cat.categories
    ncol = ncol or len(categories)
    nrow = nrow or int(np.ceil(len(categories) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow))
    axs = np.atleast_1d(axs).flatten()
    for i, cat in enumerate(categories):
        sc.pl.umap(
            adata[adata.obs[split_by] == cat],
            ax=axs[i], show=False, title=cat, **kwargs,
        )
    for j in range(i + 1, len(axs)):
        axs[j].set_visible(False)
    plt.tight_layout()


def recalculate_qc_metrics(adata):
    log.info("    Recalculating QC metrics on subset …")
    layer = "raw_counts" if "raw_counts" in adata.layers else None
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    qc_cols_to_drop = [c for c in adata.obs.columns if c in [
        "n_genes_by_counts", "total_counts", "total_counts_mt", "total_counts_hb",
        "pct_counts_mt", "pct_counts_hb",
        "pct_counts_in_top_50_genes", "pct_counts_in_top_100_genes",
        "pct_counts_in_top_200_genes", "pct_counts_in_top_500_genes",
    ]]
    adata.obs.drop(columns=qc_cols_to_drop, inplace=True, errors="ignore")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "hb"], inplace=True, log1p=False, layer=layer,
    )


def _make_barplot(crosstab_df, title, ylabel, out_path, normalize=False,
                  color_map=None):
    """Stacked barplot (black edge, legend order == stack order, config colours)."""
    save_stacked_composition(
        crosstab_df, out_path, color_map, normalize=normalize,
        ylabel=ylabel, title=title, legend_title=crosstab_df.columns.name,
        figsize=(max(8, len(crosstab_df) * 0.8), 5),
    )


def generate_barplots(adata, leiden_key, res, barplots_dir, sample_col, annot_col):
    """Generate all stacked barplot combinations for one resolution."""
    res_dir = os.path.join(barplots_dir, f"res{res}")
    os.makedirs(res_dir, exist_ok=True)

    cluster_col = leiden_key
    has_regions = (
        "region_annotation" in adata.obs.columns
        and adata.obs["region_annotation"].nunique() > 1
        and not all(adata.obs["region_annotation"] == "Unlabeled")
    )

    # 1. Cluster × Sample
    if sample_col and sample_col in adata.obs.columns:
        ct = pd.crosstab(adata.obs[cluster_col], adata.obs[sample_col])
        _make_barplot(ct, f"Cluster × Sample (res {res})", "Cells",
                      os.path.join(res_dir, "cluster_by_sample_absolute.png"))
        _make_barplot(ct, f"Cluster × Sample (res {res}) — relative", "Fraction",
                      os.path.join(res_dir, "cluster_by_sample_relative.png"), normalize=True)

        # Transposed: Sample × Cluster
        ct_t = pd.crosstab(adata.obs[sample_col], adata.obs[cluster_col])
        _make_barplot(ct_t, f"Sample × Cluster (res {res})", "Cells",
                      os.path.join(res_dir, "sample_by_cluster_absolute.png"))
        _make_barplot(ct_t, f"Sample × Cluster (res {res}) — relative", "Fraction",
                      os.path.join(res_dir, "sample_by_cluster_relative.png"), normalize=True)

    # 2. Cluster × Region
    if has_regions:
        ct = pd.crosstab(adata.obs[cluster_col], adata.obs["region_annotation"])
        _make_barplot(ct, f"Cluster × Region (res {res})", "Cells",
                      os.path.join(res_dir, "cluster_by_region_absolute.png"),
                      color_map=REGION_COLORS)
        _make_barplot(ct, f"Cluster × Region (res {res}) — relative", "Fraction",
                      os.path.join(res_dir, "cluster_by_region_relative.png"),
                      normalize=True, color_map=REGION_COLORS)

        # Transposed: Region × Cluster
        ct_t = pd.crosstab(adata.obs["region_annotation"], adata.obs[cluster_col])
        _make_barplot(ct_t, f"Region × Cluster (res {res})", "Cells",
                      os.path.join(res_dir, "region_by_cluster_absolute.png"))
        _make_barplot(ct_t, f"Region × Cluster (res {res}) — relative", "Fraction",
                      os.path.join(res_dir, "region_by_cluster_relative.png"), normalize=True)

    # 3. Sample × Region (resolution-independent, but placed per-res for context)
    if sample_col and has_regions:
        ct = pd.crosstab(adata.obs[sample_col], adata.obs["region_annotation"])
        _make_barplot(ct, f"Sample × Region", "Cells",
                      os.path.join(res_dir, "sample_by_region_absolute.png"))
        _make_barplot(ct, f"Sample × Region — relative", "Fraction",
                      os.path.join(res_dir, "sample_by_region_relative.png"), normalize=True)

    # 4. Cluster × Original annotation
    if annot_col in adata.obs.columns and adata.obs[annot_col].nunique() > 1:
        ct = pd.crosstab(adata.obs[cluster_col], adata.obs[annot_col])
        _annot_cmap = (ANNOTATION_COLORS.get(annot_col, {})
                       if isinstance(ANNOTATION_COLORS, dict) else {})
        _make_barplot(ct, f"Cluster × Annotation (res {res})", "Cells",
                      os.path.join(res_dir, "cluster_by_annotation_absolute.png"),
                      color_map=_annot_cmap)
        _make_barplot(ct, f"Cluster × Annotation (res {res}) — relative", "Fraction",
                      os.path.join(res_dir, "cluster_by_annotation_relative.png"),
                      normalize=True, color_map=_annot_cmap)


def generate_cluster_markers(adata, leiden_key, res, markers_dir, de_n_genes):
    """Generate top-N DE marker plots and TSV for one resolution."""
    res_dir = os.path.join(markers_dir, f"res{res}")
    os.makedirs(res_dir, exist_ok=True)

    n_clusters = adata.obs[leiden_key].nunique()
    if n_clusters < 2:
        log.warning("      Only %d cluster(s) at res %s. Skipping markers.", n_clusters, res)
        return

    try:
        sc.tl.rank_genes_groups(adata, groupby=leiden_key, method="wilcoxon")
        sc.tl.dendrogram(adata, groupby=leiden_key, use_rep="X_pca")
    except Exception as e:
        log.warning("      DE/dendrogram failed at res %s: %s", res, e)
        return

    # Save full results TSV
    try:
        df = sc.get.rank_genes_groups_df(adata, None)
        df.to_csv(os.path.join(res_dir, f"cluster_markers_res{res}.tsv"),
                  sep="\t", index=False)
    except Exception as e:
        log.warning("      Marker TSV failed: %s", e)

    # Dotplot
    try:
        sc.pl.rank_genes_groups_dotplot(
            adata, groupby=leiden_key, standard_scale="var",
            n_genes=de_n_genes, swap_axes=True, dendrogram=True,
        )
        plt.savefig(os.path.join(res_dir, f"dotplot_top{de_n_genes}_res{res}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        log.warning("      Dotplot failed: %s", e)

    # Heatmap
    try:
        sc.pl.rank_genes_groups_heatmap(
            adata, n_genes=de_n_genes, groupby=leiden_key, use_raw=False,
            swap_axes=True, dendrogram=True, show_gene_labels=True,
        )
        plt.savefig(os.path.join(res_dir, f"heatmap_top{de_n_genes}_res{res}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        log.warning("      Heatmap failed: %s", e)

    # Matrixplot
    try:
        sc.pl.rank_genes_groups_matrixplot(
            adata, n_genes=de_n_genes, groupby=leiden_key, use_raw=False,
            swap_axes=True, dendrogram=True,
        )
        plt.savefig(os.path.join(res_dir, f"matrixplot_top{de_n_genes}_res{res}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        log.warning("      Matrixplot failed: %s", e)


def run_clustering_branch(adata, branch_name, branch_dir, resolutions,
                          n_neighbors, subcompartment, annot_col, de_n_genes, RANDOM_SEED,
                          annotation_colors=None):
    """
    Run Leiden at multiple resolutions, produce all analysis outputs.
    """
    clust_eval_dir = os.path.join(branch_dir, "Clustering_evaluation")
    qc_dir         = os.path.join(branch_dir, "QC_plots")
    umap_dir       = os.path.join(branch_dir, "UMAPs")
    markers_dir    = os.path.join(branch_dir, "ClusterMarkers")
    barplots_dir   = os.path.join(branch_dir, "Barplots")
    for d in [clust_eval_dir, qc_dir, umap_dir, markers_dir, barplots_dir]:
        os.makedirs(d, exist_ok=True)

    # Find sample column
    sample_col = None
    for c in ["sample", "sample_batch"]:
        if c in adata.obs.columns:
            sample_col = c
            break

    # ── Leiden at all resolutions ────────────────────────────────────────
    leiden_keys = []
    for res in resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=RANDOM_SEED)
        leiden_keys.append(key)
        log.info("    [%s] res=%.1f → %d clusters", branch_name, res, adata.obs[key].nunique())

    # ── Silhouette per resolution ────────────────────────────────────────
    log.info("    [%s] Silhouette analysis …", branch_name)
    embedding_scaled = StandardScaler().fit_transform(adata.obsm["X_pca"])

    for res, key in zip(resolutions, leiden_keys):
        try:
            labels = adata.obs[key].astype(int)
            sil_values = silhouette_samples(embedding_scaled, labels)
            sil_avg = silhouette_score(embedding_scaled, labels)

            plt.figure(figsize=(12, 6))
            unique_labels = np.unique(labels)
            colours = plt.cm.tab10(np.linspace(0, 1, max(len(unique_labels), 1)))
            x_pos = 0
            for label, colour in zip(unique_labels, colours):
                vals = np.sort(sil_values[labels == label])
                xs = np.arange(len(vals)) + x_pos
                plt.bar(xs, vals, color=colour, label=f"Cluster {label}", width=1, alpha=0.7)
                x_pos = xs[-1] + 5
            plt.axhline(y=sil_avg, color="red", linestyle="--", label="Average")
            plt.ylabel("Silhouette Score")
            plt.xlabel("Cluster and Cell Index")
            plt.title(f"Silhouette — {branch_name} (res {res}) – avg {sil_avg:.3f}")
            plt.ylim(-1, 1)
            plt.legend(loc="upper right", bbox_to_anchor=(1.3, 1))
            plt.grid(axis="y")
            plt.savefig(os.path.join(clust_eval_dir, f"silhouette_res{res}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            log.warning("    [%s] Silhouette res=%.1f failed: %s", branch_name, res, e)

    # ── Clustree ─────────────────────────────────────────────────────────
    if len(leiden_keys) >= 2:
        available = [k for k in leiden_keys if k in adata.obs.columns]
        if len(available) >= 2:
            try:
                fig = clustree(adata, available, edge_weight_threshold=0.00,
                               show_fraction=True)
                fig.savefig(os.path.join(clust_eval_dir, "clustree.png"),
                            dpi=300, bbox_inches="tight")
                plt.close()
            except Exception as e:
                log.warning("    [%s] Clustree failed: %s", branch_name, e)

    # ── QC plots ─────────────────────────────────────────────────────────
    log.info("    [%s] QC plots …", branch_name)

    qc_cols = [c for c in ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"]
               if c in adata.obs.columns]
    if qc_cols:
        sc.pl.umap(adata, color=qc_cols, size=2, wspace=0.25, frameon=False)
        plt.savefig(os.path.join(qc_dir, "QC_UMAPs.png"), dpi=300, bbox_inches="tight")
        plt.close()

    violin_cols = [c for c in ["n_genes_by_counts", "total_counts",
                                "pct_counts_in_top_50_genes", "pct_counts_in_top_100_genes",
                                "pct_counts_mt", "pct_counts_hb"]
                   if c in adata.obs.columns]

    for res, key in zip(resolutions, leiden_keys):
        if violin_cols:
            try:
                layer = "raw_counts" if "raw_counts" in adata.layers else None
                sc.pl.violin(
                    adata, violin_cols, groupby=key, layer=layer,
                    jitter=0.1, multi_panel=True, show=True,
                )
                plt.savefig(os.path.join(qc_dir, f"QC_by_cluster_res{res}.png"),
                            dpi=300, bbox_inches="tight")
                plt.close()
            except Exception as e:
                log.warning("    [%s] QC violin res=%.1f failed: %s", branch_name, res, e)

    # ── UMAPs ────────────────────────────────────────────────────────────
    log.info("    [%s] UMAPs …", branch_name)

    # Apply custom palettes if configured
    if isinstance(annotation_colors, dict):
        for key in leiden_keys:
            cd = annotation_colors.get(key, {})
            if cd and key in adata.obs.columns:
                cats = adata.obs[key].cat.categories
                adata.uns[f"{key}_colors"] = [cd.get(str(c), "#cccccc") for c in cats]
        if annot_col in adata.obs.columns:
            cd = annotation_colors.get(annot_col, {})
            if cd:
                cats = adata.obs[annot_col].cat.categories
                adata.uns[f"{annot_col}_colors"] = [cd.get(str(c), "#cccccc") for c in cats]

    for res, key in zip(resolutions, leiden_keys):
        sc.pl.umap(adata, color=[key], size=2, wspace=0.25, frameon=False,
                   title=f"Leiden {res}")
        plt.savefig(os.path.join(umap_dir, f"UMAP_clusters_res{res}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    if sample_col:
        sc.pl.umap(adata, color=[sample_col], size=2, wspace=0.25, frameon=False,
                   title="By sample")
        plt.savefig(os.path.join(umap_dir, "UMAP_by_sample.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    if annot_col in adata.obs.columns:
        sc.pl.umap(adata, color=[annot_col], size=2, wspace=0.25, frameon=False,
                   title=f"Original annotation ({annot_col})")
        plt.savefig(os.path.join(umap_dir, f"UMAP_by_{annot_col}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    has_regions = (
        "region_annotation" in adata.obs.columns
        and adata.obs["region_annotation"].nunique() > 1
        and not all(adata.obs["region_annotation"] == "Unlabeled")
    )
    if has_regions:
        adata.obs["region_annotation"] = adata.obs["region_annotation"].astype("category")
        sc.pl.umap(adata, color=["region_annotation"], size=2, wspace=0.25, frameon=False,
                   title="By region")
        plt.savefig(os.path.join(umap_dir, "UMAP_by_region.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    # Split UMAPs (use middle resolution)
    log.info("    [%s] Split UMAPs …", branch_name)
    mid_idx = len(leiden_keys) // 2
    split_leiden_key = leiden_keys[mid_idx]

    if sample_col and adata.obs[sample_col].nunique() > 1:
        try:
            adata.obs[sample_col] = adata.obs[sample_col].astype("category")
            split_umap(adata, color=split_leiden_key, split_by=sample_col,
                       size=10, ncol=4)
            plt.savefig(os.path.join(umap_dir, "UMAP_split_by_sample.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            log.warning("    [%s] Split UMAP by sample failed: %s", branch_name, e)

    if has_regions:
        try:
            split_umap(adata, color=split_leiden_key, split_by="region_annotation",
                       size=10, ncol=3)
            plt.savefig(os.path.join(umap_dir, "UMAP_split_by_region.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            log.warning("    [%s] Split UMAP by region failed: %s", branch_name, e)

    # ── Cluster markers (per resolution) ─────────────────────────────────
    log.info("    [%s] Cluster markers …", branch_name)
    for res, key in zip(resolutions, leiden_keys):
        log.info("      [%s] Markers at res %.1f …", branch_name, res)
        generate_cluster_markers(adata, key, res, markers_dir, de_n_genes)

    # ── Barplots (per resolution) ────────────────────────────────────────
    log.info("    [%s] Barplots …", branch_name)
    for res, key in zip(resolutions, leiden_keys):
        generate_barplots(adata, key, res, barplots_dir, sample_col, annot_col)

    # ── Cell counts summary ──────────────────────────────────────────────
    mid_key = leiden_keys[mid_idx]
    counts_data = []
    for _, row in adata.obs.iterrows():
        entry = {"cluster": row[mid_key]}
        if sample_col:
            entry["sample"] = row[sample_col]
        if has_regions:
            entry["region"] = row["region_annotation"]
        if annot_col in adata.obs.columns:
            entry["original_annotation"] = row[annot_col]
        counts_data.append(entry)

    counts_df = pd.DataFrame(counts_data)
    summary = counts_df.groupby(list(counts_df.columns)).size().reset_index(name="n_cells")
    summary.to_csv(os.path.join(branch_dir, "cell_counts_summary.tsv"),
                   sep="\t", index=False)

    return adata


# ── Parameters ───────────────────────────────────────────────────────────────
subcompartment = snakemake.params.subcompartment
match_strings  = list(snakemake.params.strings)
annot_col      = str(snakemake.params.annot_col)
RES_MIN        = float(snakemake.params.resolution_min)
RES_MAX        = float(snakemake.params.resolution_max)
RES_STEP       = float(snakemake.params.resolution_step)
N_NEIGHBORS    = int(snakemake.params.n_neighbors)
N_PCS          = int(snakemake.params.n_pcs)
DE_N_GENES     = int(snakemake.params.de_n_genes)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS  = snakemake.params.region_colors
sub_dir        = str(snakemake.output.sub_dir)

try:
    log.info("=" * 70)
    log.info("Subclustering: %s", subcompartment)
    log.info("  Match strings: %s", match_strings)
    log.info("  Annotation col: %s", annot_col)
    log.info("  Resolutions: %.1f → %.1f (step %.1f)", RES_MIN, RES_MAX, RES_STEP)
    log.info("=" * 70)

    resolutions = _resolution_range(RES_MIN, RES_MAX, RES_STEP)
    log.info("  Resolution values: %s", resolutions)

    noharmony_dir = os.path.join(sub_dir, "NoHarmony")
    harmony_dir   = os.path.join(sub_dir, "Harmony")
    for d in [sub_dir, noharmony_dir, harmony_dir]:
        os.makedirs(d, exist_ok=True)

    # ── Load and subset ──────────────────────────────────────────────────
    adata_full = sc.read_h5ad(str(snakemake.input.adata))
    log.info("Loaded: %d cells, %d genes", adata_full.n_obs, adata_full.n_vars)

    if annot_col not in adata_full.obs.columns:
        raise ValueError(f"Annotation column '{annot_col}' not found. "
                         f"Available: {list(adata_full.obs.columns)}")

    mask = adata_full.obs[annot_col].astype(str).isin(match_strings)
    n_matched = mask.sum()
    log.info("  Matched %d / %d cells for '%s'", n_matched, adata_full.n_obs, subcompartment)

    if n_matched < 50:
        raise ValueError(f"Too few cells matched ({n_matched}). "
                         f"Check strings and annotation column.")

    adata = adata_full[mask].copy()
    del adata_full

    adata.obs["original_annotation"] = adata.obs[annot_col].copy()
    log.info("  Subset: %d cells, %d genes", adata.n_obs, adata.n_vars)

    recalculate_qc_metrics(adata)

    # ── PCA ──────────────────────────────────────────────────────────────
    log.info("  PCA …")
    sc.pp.pca(adata, use_highly_variable=False, random_state=RANDOM_SEED)

    # ══════════════════════════════════════════════════════════════════════
    # Branch 1: NoHarmony
    # ══════════════════════════════════════════════════════════════════════
    log.info("  === NoHarmony branch ===")
    adata_noharmony = adata.copy()
    sc.pp.neighbors(adata_noharmony, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_SEED)
    sc.tl.umap(adata_noharmony, random_state=RANDOM_SEED)

    adata_noharmony = run_clustering_branch(
        adata_noharmony, "NoHarmony", noharmony_dir, resolutions,
        N_NEIGHBORS, subcompartment, annot_col, DE_N_GENES, RANDOM_SEED,
        annotation_colors=ANNOTATION_COLORS,
    )

    log.info("  Saving NoHarmony adata …")
    adata_noharmony.write(os.path.join(noharmony_dir,
                                        f"adata_{subcompartment}_noharmony.h5ad"))

    # ══════════════════════════════════════════════════════════════════════
    # Branch 2: Harmony
    # ══════════════════════════════════════════════════════════════════════
    log.info("  === Harmony branch ===")
    adata_harmony = adata.copy()

    sample_col = None
    for c in ["sample", "sample_batch"]:
        if c in adata_harmony.obs.columns:
            sample_col = c
            break

    if sample_col and adata_harmony.obs[sample_col].nunique() > 1:
        adata_harmony.obsm["X_pca_original"] = adata_harmony.obsm["X_pca"].copy()
        sc.external.pp.harmony_integrate(adata_harmony, key=sample_col, random_state=RANDOM_SEED)
        adata_harmony.obsm["X_pca"] = adata_harmony.obsm["X_pca_harmony"]
    else:
        log.warning("  Only 1 sample or no sample column — Harmony skipped.")

    sc.pp.neighbors(adata_harmony, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_SEED)
    sc.tl.umap(adata_harmony, random_state=RANDOM_SEED)

    adata_harmony = run_clustering_branch(
        adata_harmony, "Harmony", harmony_dir, resolutions,
        N_NEIGHBORS, subcompartment, annot_col, DE_N_GENES, RANDOM_SEED,
        annotation_colors=ANNOTATION_COLORS,
    )

    log.info("  Saving Harmony adata …")
    adata_harmony.write(os.path.join(harmony_dir,
                                      f"adata_{subcompartment}_harmony.h5ad"))

    log.info("Subclustering complete for %s.", subcompartment)

except Exception:
    log.error("FAILED for %s:\n%s", subcompartment, traceback.format_exc())
    raise
