"""
leiden_analysis.py – Per-resolution visualisation (plots only)
===============================================================
Reads the preprocessed adata which already has leiden_{res} in obs
(computed in preprocess_umap).  Produces all plots and the marker
TSV.  No h5ad is saved — avoids duplicating the large object.
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
from pyclustree import clustree
import scanpy as sc
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.preprocessing import StandardScaler


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
log = logging.getLogger("leiden_analysis")

sample_id      = snakemake.params.sample_id
resolution     = float(snakemake.params.resolution)
DE_N_GENES     = int(snakemake.params.de_n_genes)
RES_SCAN_MIN   = float(snakemake.params.resolution_scan_min)
RES_SCAN_MAX   = float(snakemake.params.resolution_scan_max)
RES_SCAN_STEP  = float(snakemake.params.resolution_scan_step)
markers_path   = str(snakemake.input.cell_markers)
res_dir        = str(snakemake.output.res_dir)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors


def read_tsv_to_dict(tsv_path):
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = {field: [] for field in reader.fieldnames}
        for row in reader:
            for key, val in row.items():
                if val.strip():
                    columns[key].append(val.strip())
    return columns


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


try:
    log.info("=" * 70)
    log.info("Leiden analysis: sample=%s, resolution=%.2f", sample_id, resolution)
    log.info("=" * 70)

    # Create subdirectories
    clust_eval_dir  = os.path.join(res_dir, "clustering_evaluation")
    qc_dir          = os.path.join(res_dir, "QC_plots")
    spatial_dir     = os.path.join(res_dir, "spatial_clusters")
    umap_dir        = os.path.join(res_dir, "UMAPs")
    markers_out_dir = os.path.join(res_dir, "Markers")
    for d in [res_dir, clust_eval_dir, qc_dir, spatial_dir, umap_dir, markers_out_dir]:
        os.makedirs(d, exist_ok=True)

    # ── Load adata (leiden columns already computed) ─────────────────────
    adata = sc.read_h5ad(str(snakemake.input.adata))

    # Use the pre-computed leiden column for this resolution
    leiden_key = f"leiden_{str(resolution).replace('.', '_')}"
    if leiden_key not in adata.obs.columns:
        log.warning("Column '%s' not found, computing on the fly …", leiden_key)
        sc.tl.leiden(adata, resolution=resolution, key_added=leiden_key)

    # Set as the working "leiden" column for plotting convenience
    adata.obs["leiden"] = adata.obs[leiden_key]
    # Ensure string categories with numeric ordering preserved (precomputed
    # clusters may be int, which breaks sc.pl.rank_genes_groups_dotplot)
    raw_cats = adata.obs["leiden"].cat.categories
    sorted_str_cats = [str(c) for c in sorted(raw_cats, key=lambda x: int(x))]
    adata.obs["leiden"] = (adata.obs["leiden"].astype(str)
                           .astype(pd.CategoricalDtype(categories=sorted_str_cats,
                                                        ordered=False)))
    log.info("Using %s: %d clusters", leiden_key, adata.obs["leiden"].nunique())

    # Apply custom palette if configured
    def _set_palette(adata, obs_key, colors_cfg):
        if not isinstance(colors_cfg, dict):
            return
        cd = colors_cfg.get(obs_key, {})
        if cd and obs_key in adata.obs.columns:
            cats = adata.obs[obs_key].cat.categories
            adata.uns[f"{obs_key}_colors"] = [cd.get(str(c), "#cccccc") for c in cats]

    _set_palette(adata, "leiden", ANNOTATION_COLORS)

    # Apply custom region palette if configured
    if REGION_COLORS and "region_annotation" in adata.obs.columns:
        adata.obs["region_annotation"] = adata.obs["region_annotation"].astype("category")
        cats = adata.obs["region_annotation"].cat.categories
        adata.uns["region_annotation_colors"] = [
            REGION_COLORS.get(str(c), "#cccccc") for c in cats
        ]

    has_annotations = (
        "region_annotation" in adata.obs.columns
        and adata.obs["region_annotation"].nunique() > 1
        and not all(adata.obs["region_annotation"] == "Unlabeled")
    )

    # ── Clustering evaluation (silhouette) ───────────────────────────────
    log.info("Silhouette analysis …")
    embedding_scaled = StandardScaler().fit_transform(adata.obsm["X_pca"])
    labels = adata.obs["leiden"].astype(int)
    sil_values = silhouette_samples(embedding_scaled, labels)
    sil_avg = silhouette_score(embedding_scaled, labels)

    plt.figure(figsize=(12, 6))
    unique_labels = np.unique(labels)
    colours = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))
    x_pos = 0
    for label, colour in zip(unique_labels, colours):
        vals = np.sort(sil_values[labels == label])
        xs = np.arange(len(vals)) + x_pos
        plt.bar(xs, vals, color=colour, label=f"Cluster {label}", width=1, alpha=0.7)
        x_pos = xs[-1] + 5
    plt.axhline(y=sil_avg, color="red", linestyle="--", label="Average Silhouette")
    plt.ylabel("Silhouette Score")
    plt.xlabel("Cluster and Cell Index")
    plt.title(f"Silhouette (resolution {resolution}) – avg {sil_avg:.3f}")
    plt.ylim(-1, 1)
    plt.legend(loc="upper right", bbox_to_anchor=(1.3, 1))
    plt.grid(axis="y")
    plt.savefig(os.path.join(clust_eval_dir, f"silhouette_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # Clustree (uses all pre-computed leiden columns)
    resolutions = np.arange(RES_SCAN_MIN, RES_SCAN_MAX, RES_SCAN_STEP).round(1)
    leiden_keys = [f"leiden_{str(r).replace('.', '_')}" for r in resolutions]
    available_keys = [k for k in leiden_keys if k in adata.obs.columns]
    if len(available_keys) >= 2:
        fig = clustree(adata, available_keys, edge_weight_threshold=0.00, show_fraction=True)
        fig.savefig(os.path.join(clust_eval_dir, f"clustree_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    # ── QC plots ─────────────────────────────────────────────────────────
    log.info("QC plots …")
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts",
         "pct_counts_in_top_50_genes", "pct_counts_in_top_100_genes",
         "pct_counts_mt", "pct_counts_hb"],
        groupby="leiden", layer="raw_counts", jitter=0.1, multi_panel=True, show=True,
    )
    plt.savefig(os.path.join(qc_dir, f"QC_by_cluster_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(
        adata,
        color=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"],
        size=2, wspace=0.25, frameon=False,
    )
    plt.savefig(os.path.join(qc_dir, f"QC_UMAPs_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    if has_annotations:
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts",
             "pct_counts_in_top_50_genes", "pct_counts_in_top_100_genes"],
            groupby="region_annotation", layer="raw_counts", jitter=0.1,
            multi_panel=True, show=True,
        )
        plt.savefig(os.path.join(qc_dir, f"QC_by_region_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    # ── UMAPs ────────────────────────────────────────────────────────────
    log.info("UMAP plots …")
    sc.pl.umap(adata, color=["leiden"], size=2, wspace=0.25, frameon=False)
    plt.savefig(os.path.join(res_dir, f"UMAP_clusters_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    if has_annotations:
        sc.pl.umap(adata, color=["region_annotation"], size=2, wspace=0.25, frameon=False)
        plt.savefig(os.path.join(umap_dir, f"UMAP_region_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

        split_umap(adata, color="leiden", split_by="region_annotation", size=10)
        plt.savefig(os.path.join(umap_dir, f"UMAP_split_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

        ax = (pd.crosstab(adata.obs["leiden"], adata.obs["region_annotation"],
                          normalize="columns").T.plot(kind="bar", stacked=True))
        ax.legend(title=f"leiden_{resolution}", bbox_to_anchor=(1.26, 1.02), loc="upper right")
        plt.savefig(os.path.join(umap_dir, f"normalized_barplot_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

        ax = (pd.crosstab(adata.obs["leiden"], adata.obs["region_annotation"])
              .T.plot(kind="bar", stacked=True))
        ax.legend(title=f"leiden_{resolution}", bbox_to_anchor=(1.26, 1.02), loc="upper right")
        plt.savefig(os.path.join(umap_dir, f"absolute_barplot_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    # ── Spatial cluster maps ─────────────────────────────────────────────
    log.info("Spatial plots …")
    sc.pl.spatial(adata, color="leiden", spot_size=20, title="Leiden Clusters", frameon=False)
    plt.savefig(os.path.join(res_dir, f"spatial_all_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    for cluster in adata.obs["leiden"].unique():
        adata_sub = adata[adata.obs["leiden"] == cluster, :]
        sc.pl.spatial(
            adata_sub, color="leiden", spot_size=20,
            title=f"Cluster {cluster}", frameon=False,
            palette=["black"], alpha_img=0.5,
        )
        plt.savefig(os.path.join(spatial_dir, f"spatial_{cluster}_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    for cluster in adata.obs["leiden"].cat.categories:
        adata.obs["_highlight"] = pd.Categorical(
            adata.obs["leiden"].apply(
                lambda x, c=cluster: str(c) if x == c else "Other"
            ),
            categories=["Other", str(cluster)],
        )
        sc.pl.umap(
            adata, color="_highlight", size=2, wspace=0.25, frameon=False,
            palette={"Other": "gray", str(cluster): "red"},
            title=f"Cluster {cluster}",
        )
        plt.savefig(os.path.join(spatial_dir, f"UMAP_highlight_{cluster}_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    # ── Markers ──────────────────────────────────────────────────────────
    log.info("DE + marker visualisation …")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    sc.tl.dendrogram(adata, groupby="leiden", use_rep="X_pca")

    sc.pl.rank_genes_groups_dotplot(
        adata, groupby="leiden", standard_scale="var", n_genes=DE_N_GENES,
        swap_axes=True, dendrogram=True,
    )
    plt.savefig(os.path.join(res_dir, f"top{DE_N_GENES}_dotplot_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.rank_genes_groups_heatmap(
        adata, n_genes=DE_N_GENES, groupby="leiden", use_raw=False,
        swap_axes=True, dendrogram=True, show_gene_labels=True,
    )
    plt.savefig(os.path.join(res_dir, f"top{DE_N_GENES}_heatmap_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.rank_genes_groups_matrixplot(
        adata, n_genes=DE_N_GENES, groupby="leiden", use_raw=False,
        swap_axes=True, dendrogram=True,
    )
    plt.savefig(os.path.join(res_dir, f"top{DE_N_GENES}_matrixplot_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # Cell-type markers
    marker_dict = read_tsv_to_dict(markers_path)
    marker_dict = {
        ct: [g for g in genes if g in adata.var_names]
        for ct, genes in marker_dict.items()
    }
    marker_dict = {ct: genes for ct, genes in marker_dict.items() if genes}

    for cell_type, genes in marker_dict.items():
        ct_dir = os.path.join(markers_out_dir, cell_type)
        os.makedirs(ct_dir, exist_ok=True)

        sc.pl.dotplot(
            adata, genes, groupby="leiden", dendrogram=True,
            standard_scale="var", show=True, title=f"{cell_type} markers",
        )
        plt.savefig(os.path.join(ct_dir, f"dotplot_{cell_type}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

        for gene in genes:
            sc.pl.umap(adata, color=gene, show=False)
            plt.savefig(os.path.join(ct_dir, f"UMAP_{gene}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()

    # ── Save cluster markers TSV only ────────────────────────────────────
    df = sc.get.rank_genes_groups_df(adata, None)
    df.to_csv(os.path.join(res_dir, f"cluster_markers_{sample_id}.tsv"), sep="\t", index=False)

    log.info("Leiden analysis complete: sample=%s, res=%.2f", sample_id, resolution)

except Exception:
    log.error("FAILED: sample=%s, res=%.2f\n%s", sample_id, resolution, traceback.format_exc())
    raise
