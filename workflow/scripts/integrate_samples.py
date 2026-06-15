"""
integrate_samples.py – Multi-sample integration
=================================================
1. Concatenate all annotated samples
2. PCA + UMAP (no batch correction) → concatenated.h5ad
3. Harmony integration + UMAP → harmony_integrated.h5ad
4. Geosketch 25% of the concatenated data, cluster the sketch,
   ingest-project labels onto the full dataset → sketched.h5ad
"""

import gc
import logging
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scanpy.external as sce

# Shared composition-barplot helpers (black edge, legend==stack order, config colours)
try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:                      # very old Snakemake
    _here = os.getcwd()
sys.path.insert(0, _here)
from composition_barplots import (
    composition_pair, find_niche_column, composition_grouped,
)

try:
    import geosketch as sketch
    HAS_GEOSKETCH = True
except ImportError:
    HAS_GEOSKETCH = False


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
log = logging.getLogger("integrate_samples")

# ── Parameters ───────────────────────────────────────────────────────────────
annotated_paths = [str(p) for p in snakemake.input.annotated]
sample_ids      = list(snakemake.params.sample_ids)
N_NEIGHBORS     = int(snakemake.params.n_neighbors)
SKETCH_FRAC     = float(snakemake.params.sketch_fraction)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors
DPI          = int(getattr(snakemake.params, "dpi", 300))
NICHE_COLUMN      = getattr(snakemake.params, "niche_column", "")

out_concat     = str(snakemake.output.concatenated)
out_harmony    = str(snakemake.output.harmony)
out_sketched   = str(snakemake.output.sketched)
output_dir     = str(Path(out_concat).parent)

try:
    log.info("=" * 70)
    log.info("Integrating %d samples", len(annotated_paths))
    log.info("=" * 70)

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # ── 1. Load and concatenate ──────────────────────────────────────────
    adatas = []
    for path, sid in zip(annotated_paths, sample_ids):
        log.info("  Loading %s …", sid)
        ad = sc.read_h5ad(path)
        # Ensure sample label is present
        if "sample" not in ad.obs.columns:
            ad.obs["sample"] = sid
        adatas.append(ad)

    log.info("Concatenating …")
    adata = sc.concat(adatas, join="inner", label="sample_batch", keys=sample_ids)
    log.info("  Concatenated: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # Carry over raw_counts if available
    if "raw_counts" in adatas[0].layers:
        # sc.concat preserves layers if they exist in all objects
        if "raw_counts" not in adata.layers:
            log.warning("  raw_counts layer lost during concatenation")

    del adatas
    gc.collect()

    # ── 2. PCA + UMAP (no batch correction) ──────────────────────────────
    log.info("PCA + UMAP (uncorrected) …")
    sc.pp.pca(adata, use_highly_variable=False)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata)

    # Leiden on uncorrected
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_uncorrected")

    log.info("Saving concatenated → %s", out_concat)
    adata.write(out_concat)

    # Apply custom palettes if configured
    if isinstance(ANNOTATION_COLORS, dict):
        for obs_key in ["sample_batch", "cell_type_tsv", "cell_type_refined",
                        "cell_type_ingest", "cell_type_external", "leiden_uncorrected"]:
            cd = ANNOTATION_COLORS.get(obs_key, {})
            if cd and obs_key in adata.obs.columns:
                cats = adata.obs[obs_key].cat.categories if hasattr(adata.obs[obs_key], "cat") else []
                if len(cats) > 0:
                    adata.uns[f"{obs_key}_colors"] = [cd.get(str(c), "#cccccc") for c in cats]

    # Plot uncorrected
    sc.pl.umap(adata, color=["sample_batch"], size=2, frameon=False,
               title="Uncorrected – by sample")
    plt.savefig(os.path.join(output_dir, "UMAP_uncorrected_by_sample.png"),
                dpi=DPI, bbox_inches="tight")
    plt.close()

    if "cell_type_tsv" in adata.obs.columns:
        sc.pl.umap(adata, color=["cell_type_tsv"], size=2, frameon=False,
                   title="Uncorrected – by cell type (TSV)")
        plt.savefig(os.path.join(output_dir, "UMAP_uncorrected_by_celltype.png"),
                    dpi=DPI, bbox_inches="tight")
        plt.close()

    # ── 3. Harmony integration ───────────────────────────────────────────
    log.info("Harmony integration …")
    adata_harmony = adata.copy()
    adata_harmony.obsm["X_pca_uncorrected"] = adata_harmony.obsm["X_pca"].copy()

    sc.external.pp.harmony_integrate(adata_harmony, key="sample_batch")
    adata_harmony.obsm["X_pca"] = adata_harmony.obsm["X_pca_harmony"]

    sc.pp.neighbors(adata_harmony, n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata_harmony)
    sc.tl.leiden(adata_harmony, resolution=0.5, key_added="leiden_harmony")

    log.info("Saving Harmony integrated → %s", out_harmony)
    adata_harmony.write(out_harmony)

    # Copy palettes to harmony object
    for key in list(adata.uns.keys()):
        if key.endswith("_colors") and key not in adata_harmony.uns:
            adata_harmony.uns[key] = adata.uns[key]

    # Plots
    sc.pl.umap(adata_harmony, color=["sample_batch"], size=2, frameon=False,
               title="Harmony – by sample")
    plt.savefig(os.path.join(output_dir, "UMAP_harmony_by_sample.png"),
                dpi=DPI, bbox_inches="tight")
    plt.close()

    if "cell_type_tsv" in adata_harmony.obs.columns:
        sc.pl.umap(adata_harmony, color=["cell_type_tsv"], size=2, frameon=False,
                   title="Harmony – by cell type (TSV)")
        plt.savefig(os.path.join(output_dir, "UMAP_harmony_by_celltype.png"),
                    dpi=DPI, bbox_inches="tight")
        plt.close()

    if "cell_type_refined" in adata_harmony.obs.columns:
        sc.pl.umap(adata_harmony, color=["cell_type_refined"], size=2, frameon=False,
                   title="Harmony – by cell type (refined)")
        plt.savefig(os.path.join(output_dir, "UMAP_harmony_by_celltype_refined.png"),
                    dpi=DPI, bbox_inches="tight")
        plt.close()

    # ── Composition barplots (sample / region / niche) ───────────────────
    log.info("Composition barplots …")
    bar_dir = os.path.join(output_dir, "barplots")
    os.makedirs(bar_dir, exist_ok=True)
    annot_cols = [c for c in ["cell_type_tsv", "cell_type_refined",
                              "cell_type_ingest", "cell_type_external"]
                  if c in adata.obs.columns]
    has_regions = ("region_annotation" in adata.obs.columns
                   and adata.obs["region_annotation"].nunique() > 1
                   and not all(adata.obs["region_annotation"] == "Unlabeled"))
    niche_col = find_niche_column(adata, NICHE_COLUMN)

    for annot_col in annot_cols:
        cmap = (ANNOTATION_COLORS.get(annot_col, {})
                if isinstance(ANNOTATION_COLORS, dict) else {})
        label = annot_col.replace("cell_type_", "")
        composition_pair(adata, "sample_batch", annot_col, cmap, bar_dir,
                         f"barplot_{label}_by_sample", group_label="Sample",
                         cat_label=label, dpi=DPI)
        if has_regions:
            composition_pair(adata, "region_annotation", annot_col, cmap,
                             bar_dir, f"barplot_{label}_by_region",
                             group_label="Region", cat_label=label, dpi=DPI)
        if niche_col:
            composition_pair(adata, niche_col, annot_col, cmap, bar_dir,
                             f"barplot_{label}_by_niche", group_label="Niche",
                             cat_label=label, dpi=DPI)

    # Sample × region cell-type composition — two single-axis grouped layouts,
    # absolute + relative:
    #   • grouped by sample (sample headers, regions within each group)
    #   • grouped by region (region headers, samples within each group)
    if has_regions and annot_cols:
        primary = annot_cols[0]
        pcmap = (ANNOTATION_COLORS.get(primary, {})
                 if isinstance(ANNOTATION_COLORS, dict) else {})
        plabel = primary.replace("cell_type_", "")
        region_order = (list(REGION_COLORS.keys())
                        if isinstance(REGION_COLORS, dict) and REGION_COLORS else None)
        sample_order = sorted(adata.obs["sample_batch"].astype(str).unique())
        for norm, suffix, yl in [(False, "absolute", "Number of cells"),
                                 (True, "relative", "Percentage (%)")]:
            composition_grouped(
                adata, "sample_batch", "region_annotation", primary, pcmap,
                os.path.join(bar_dir, f"barplot_{plabel}_grouped_by_sample_{suffix}.png"),
                normalize=norm, outer_order=sample_order, inner_order=region_order,
                cat_label=plabel, ylabel=yl, dpi=DPI)
            composition_grouped(
                adata, "region_annotation", "sample_batch", primary, pcmap,
                os.path.join(bar_dir, f"barplot_{plabel}_grouped_by_region_{suffix}.png"),
                normalize=norm, outer_order=region_order, inner_order=sample_order,
                cat_label=plabel, ylabel=yl, dpi=DPI)

    # ── 4. Geosketch ─────────────────────────────────────────────────────
    if not HAS_GEOSKETCH:
        log.warning("geosketch not installed — skipping sketching. "
                     "Saving Harmony object as sketch placeholder.")
        adata_harmony.write(out_sketched)
    else:
        log.info("Geosketching %.0f%% of %d cells …", SKETCH_FRAC * 100, adata.n_obs)
        n_sketch = max(100, int(adata.n_obs * SKETCH_FRAC))

        # Sketch from PCA embedding
        pca_data = adata.obsm["X_pca"]
        sketch_indices = sketch.gs(pca_data, n_sketch, replace=False)

        sketched_adata = adata[sketch_indices].copy()
        log.info("  Sketched: %d cells", sketched_adata.n_obs)

        # Cluster the sketch
        sc.pp.neighbors(sketched_adata, n_neighbors=N_NEIGHBORS)
        sc.tl.umap(sketched_adata)
        sc.tl.leiden(sketched_adata, resolution=0.5, key_added="clusters")

        # Project labels back to full dataset via ingest (in chunks)
        log.info("  Projecting sketch clusters onto full dataset …")
        chunk_size = 100000
        n_chunks = (adata.n_obs + chunk_size - 1) // chunk_size
        ingested_chunks = []

        for i in range(n_chunks):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, adata.n_obs)
            chunk = adata[start:end].copy()
            sc.tl.ingest(chunk, sketched_adata, obs="clusters")
            ingested_chunks.append(chunk)
            log.info("    Ingested chunk %d/%d (%d cells)", i + 1, n_chunks, chunk.n_obs)

        adata_ingested = sc.concat(ingested_chunks)

        log.info("Saving sketched → %s", out_sketched)
        adata_ingested.write(out_sketched)

        # Plot
        sc.pl.umap(adata_ingested, color=["clusters"], size=2, frameon=False,
                   title="Geosketch – ingested clusters")
        plt.savefig(os.path.join(output_dir, "UMAP_sketched_clusters.png"),
                    dpi=DPI, bbox_inches="tight")
        plt.close()

        sc.pl.umap(adata_ingested, color=["sample_batch"], size=2, frameon=False,
                   title="Geosketch – by sample")
        plt.savefig(os.path.join(output_dir, "UMAP_sketched_by_sample.png"),
                    dpi=DPI, bbox_inches="tight")
        plt.close()

        del ingested_chunks, sketched_adata
        gc.collect()

    log.info("Integration complete.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
