"""
ingest_ref.py – Transfer cell-type labels from a reference via sc.tl.ingest
=============================================================================
Reads the preprocessed query adata and a reference h5ad with cell-type
annotations.  Projects the query onto the reference PCA space and
transfers labels using k-NN classification (scanpy ingest).

Produces:
- obs['cell_type_ingest']: transferred labels
- UMAP and spatial plots of the ingested annotations
- Top-N DE markers for the ingested cell types
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
log = logging.getLogger("ingest_ref")

sample_id      = snakemake.params.sample_id
ref_label_key  = snakemake.params.ref_label_key
DE_N_GENES     = int(snakemake.params.de_n_genes)
ANNOTATION_COLORS = snakemake.params.annotation_colors

adata_path     = str(snakemake.input.adata)
ref_path       = str(snakemake.input.ingest_ref)
out_adata_path = str(snakemake.output.adata_ingested)
plots_dir      = str(snakemake.output.plots_dir)

try:
    log.info("=" * 70)
    log.info("Ingest reference: sample=%s", sample_id)
    log.info("  Reference: %s", ref_path)
    log.info("  Label key: %s", ref_label_key)
    log.info("=" * 70)

    os.makedirs(plots_dir, exist_ok=True)

    # ── Load query and reference ─────────────────────────────────────────
    log.info("Loading query adata …")
    adata = sc.read_h5ad(adata_path)
    log.info("  Query: %d cells, %d genes", adata.n_obs, adata.n_vars)

    log.info("Loading reference adata …")
    adata_ref = sc.read_h5ad(ref_path)
    log.info("  Reference: %d cells, %d genes", adata_ref.n_obs, adata_ref.n_vars)

    # Validate reference
    if ref_label_key not in adata_ref.obs.columns:
        available = [c for c in adata_ref.obs.columns
                     if adata_ref.obs[c].dtype == "category" or adata_ref.obs[c].dtype == "object"]
        raise ValueError(
            f"Label key '{ref_label_key}' not found in reference. "
            f"Available categorical columns: {available}"
        )

    if "X_pca" not in adata_ref.obsm:
        log.info("  Reference has no PCA — computing …")
        sc.pp.pca(adata_ref)

    adata_ref.obs[ref_label_key] = adata_ref.obs[ref_label_key].astype("category")
    n_ref_types = adata_ref.obs[ref_label_key].nunique()
    log.info("  Reference has %d cell types: %s",
             n_ref_types, list(adata_ref.obs[ref_label_key].cat.categories))

    # ── Subset to shared genes ───────────────────────────────────────────
    shared_genes = adata.var_names.intersection(adata_ref.var_names)
    log.info("  Shared genes: %d (query: %d, ref: %d)",
             len(shared_genes), adata.n_vars, adata_ref.n_vars)

    if len(shared_genes) < 100:
        raise ValueError(
            f"Only {len(shared_genes)} shared genes — too few for ingest. "
            f"Check that query and reference use the same gene annotation."
        )

    adata_ref = adata_ref[:, shared_genes].copy()
    adata_query = adata[:, shared_genes].copy()

    # Recompute PCA on shared genes for the reference
    log.info("  Recomputing reference PCA on shared genes …")
    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)

    # ── Ingest ───────────────────────────────────────────────────────────
    log.info("Running sc.tl.ingest …")
    sc.tl.ingest(adata_query, adata_ref, obs=ref_label_key)
    log.info("  Ingest complete.")

    # Transfer the label back to the full adata (all genes)
    adata.obs["cell_type_ingest"] = adata_query.obs[ref_label_key].values
    adata.obs["cell_type_ingest"] = adata.obs["cell_type_ingest"].astype("category")

    # Also store the ingest UMAP (computed on shared genes)
    adata.obsm["X_umap_ingest"] = adata_query.obsm["X_umap"]

    ingest_counts = adata.obs["cell_type_ingest"].value_counts()
    log.info("Ingest annotation distribution:\n%s", ingest_counts.to_string())

    # ── Plots ────────────────────────────────────────────────────────────
    log.info("Generating ingest plots …")

    # Apply custom palette if configured
    if isinstance(ANNOTATION_COLORS, dict):
        cd = ANNOTATION_COLORS.get("cell_type_ingest", {})
        if cd:
            cats = adata.obs["cell_type_ingest"].cat.categories
            adata.uns["cell_type_ingest_colors"] = [cd.get(str(c), "#cccccc") for c in cats]

    library_id = None
    if "spatial" in adata.uns and len(adata.uns["spatial"]) == 1:
        library_id = list(adata.uns["spatial"].keys())[0]

    # UMAP (original embedding) coloured by ingest labels
    sc.pl.umap(adata, color=["cell_type_ingest"], size=2, frameon=False,
               title="Ingest annotation (original UMAP)")
    plt.savefig(os.path.join(plots_dir, f"UMAP_ingest_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # UMAP from ingest projection
    adata.obsm["X_umap_orig"] = adata.obsm["X_umap"].copy()
    adata.obsm["X_umap"] = adata.obsm["X_umap_ingest"]
    sc.pl.umap(adata, color=["cell_type_ingest"], size=2, frameon=False,
               title="Ingest annotation (projected UMAP)")
    plt.savefig(os.path.join(plots_dir, f"UMAP_ingest_projected_{sample_id}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    # Restore original UMAP
    adata.obsm["X_umap"] = adata.obsm["X_umap_orig"]
    del adata.obsm["X_umap_orig"]

    # Spatial overview
    try:
        sc.pl.spatial(adata, color="cell_type_ingest", spot_size=20,
                      frameon=False, title="Ingest – spatial",
                      library_id=library_id)
        plt.savefig(os.path.join(plots_dir, f"spatial_ingest_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        log.warning("Spatial plot failed: %s", e)

    # Per-cell-type spatial
    spatial_ct_dir = os.path.join(plots_dir, "spatial_per_celltype")
    os.makedirs(spatial_ct_dir, exist_ok=True)
    for ct in adata.obs["cell_type_ingest"].cat.categories:
        adata.obs["_hl"] = adata.obs["cell_type_ingest"].apply(
            lambda x, c=ct: c if x == c else "Other"
        )
        try:
            sc.pl.spatial(adata, color="_hl", spot_size=20, frameon=False,
                          palette={"Other": "#d3d3d3", ct: "#e41a1c"},
                          title=ct, library_id=library_id)
            safe = ct.replace("/", "_").replace(" ", "_")
            plt.savefig(os.path.join(spatial_ct_dir, f"{safe}_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            log.warning("Spatial '%s' failed: %s", ct, e)
    if "_hl" in adata.obs.columns:
        del adata.obs["_hl"]

    # DE markers
    valid_types = ingest_counts[ingest_counts >= 10]
    if len(valid_types) >= 2:
        log.info("DE markers (top %d) …", DE_N_GENES)
        sub = adata[adata.obs["cell_type_ingest"].isin(valid_types.index)].copy()
        sub.obs["cell_type_ingest"] = sub.obs["cell_type_ingest"].cat.remove_unused_categories()
        try:
            sc.tl.rank_genes_groups(sub, groupby="cell_type_ingest", method="wilcoxon")
            sc.tl.dendrogram(sub, groupby="cell_type_ingest", use_rep="X_pca")

            sc.pl.rank_genes_groups_dotplot(
                sub, groupby="cell_type_ingest", standard_scale="var",
                n_genes=DE_N_GENES, swap_axes=True, dendrogram=True)
            plt.savefig(os.path.join(plots_dir, f"top{DE_N_GENES}_dotplot_ingest_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()

            sc.pl.rank_genes_groups_heatmap(
                sub, n_genes=DE_N_GENES, groupby="cell_type_ingest", use_raw=False,
                swap_axes=True, dendrogram=True, show_gene_labels=True)
            plt.savefig(os.path.join(plots_dir, f"top{DE_N_GENES}_heatmap_ingest_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()

            sc.pl.rank_genes_groups_matrixplot(
                sub, n_genes=DE_N_GENES, groupby="cell_type_ingest", use_raw=False,
                swap_axes=True, dendrogram=True)
            plt.savefig(os.path.join(plots_dir, f"top{DE_N_GENES}_matrixplot_ingest_{sample_id}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close()

            df = sc.get.rank_genes_groups_df(sub, None)
            df.to_csv(os.path.join(plots_dir, f"DE_markers_ingest_{sample_id}.tsv"),
                      sep="\t", index=False)
        except Exception as e:
            log.warning("DE failed: %s", e)

    # ── Save ─────────────────────────────────────────────────────────────
    log.info("Saving ingested adata → %s", out_adata_path)
    Path(out_adata_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_adata_path)

    log.info("Ingest complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
