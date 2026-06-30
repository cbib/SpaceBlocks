"""
preprocess_umap.py – CORE: QC, normalise, embed, cluster.
=========================================================
Reads the UNFILTERED contract h5ad (raw counts in X, obsm["spatial"],
obs["region_annotation"], optional uns["spatial"] image) produced by the
technology-specific head (prepare_input) and validated upstream. Technology-
specific I/O (Space Ranger, image building, geojson) lives in the head; this rule
is the common core: QC filter → normalise → PCA → UMAP → multi-resolution Leiden,
with optional reload of precomputed clusters. Saves the per-sample h5ad and a
metadata TSV.
"""
import logging
import os
import sys
import traceback
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S", handlers=log_handlers)
log = logging.getLogger("preprocess_umap")

# ── Parameters ───────────────────────────────────────────────────────────────
in_h5ad        = str(snakemake.input.h5ad)
sample_id      = snakemake.params.sample_id
MIN_COUNTS     = int(snakemake.params.min_counts)
MIN_CELLS      = int(snakemake.params.min_cells)
MIN_GENES      = int(snakemake.params.min_genes)
N_NEIGHBORS    = int(snakemake.params.n_neighbors)
N_PCS          = int(snakemake.params.n_pcs)
RES_SCAN_MIN   = float(snakemake.params.resolution_scan_min)
RES_SCAN_MAX   = float(snakemake.params.resolution_scan_max)
RES_SCAN_STEP  = float(snakemake.params.resolution_scan_step)
RANDOM_SEED    = int(snakemake.params.random_seed)
USE_PRECOMPUTED = bool(snakemake.params.use_precomputed)
PRECOMPUTED_DIR = str(snakemake.params.precomputed_metadata_dir)
REGION_COLORS   = snakemake.params.region_colors

out_adata      = str(snakemake.output.adata)
out_metadata   = str(snakemake.output.metadata)
output_dir     = str(Path(out_adata).parent)


def _resolution_range(rmin, rmax, step):
    n = round((rmax - rmin) / step)
    return [round(rmin + i * step, 10) for i in range(n + 1)]


try:
    log.info("=" * 70)
    log.info("Preprocessing (core) sample: %s", sample_id)
    log.info("  input: %s", in_h5ad)
    log.info("=" * 70)

    # ── 1. Load the validated, unfiltered contract h5ad ──────────────────
    adata = sc.read_h5ad(in_h5ad)
    if "sample" not in adata.obs.columns:
        adata.obs["sample"] = sample_id
    log.info("  loaded %d cells, %d genes (raw, unfiltered)", adata.n_obs, adata.n_vars)

    # ── 2. QC metrics + filter ───────────────────────────────────────────
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hb"], inplace=True, log1p=False)

    log.info("QC: min_counts=%d, min_cells=%d, min_genes=%d", MIN_COUNTS, MIN_CELLS, MIN_GENES)
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    log.info("After QC: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # ── 3. Region annotation plot (already joined upstream) ──────────────
    if "region_annotation" not in adata.obs.columns:
        adata.obs["region_annotation"] = "Unlabeled"
    if isinstance(REGION_COLORS, dict) and REGION_COLORS:
        adata.obs["region_annotation"] = adata.obs["region_annotation"].astype("category")
        cats = adata.obs["region_annotation"].cat.categories
        adata.uns["region_annotation_colors"] = [REGION_COLORS.get(str(c), "#cccccc") for c in cats]
    try:
        spatial_kw = ({"library_id": list(adata.uns["spatial"].keys())[0]}
                      if isinstance(adata.uns.get("spatial"), dict) and adata.uns["spatial"]
                      else {"spot_size": 20})
        sc.pl.spatial(adata, color="region_annotation", title="Annotated Regions",
                      show=False, **spatial_kw)
        plt.savefig(os.path.join(output_dir, f"region_annotation_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        log.warning("Region plot failed: %s", e)

    # ── 4. Normalise ─────────────────────────────────────────────────────
    log.info("Normalising (target_sum=None + log1p) …")
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)

    # ── 5. PCA → UMAP ────────────────────────────────────────────────────
    log.info("PCA (all genes) → UMAP (seed=%d) …", RANDOM_SEED)
    sc.pp.pca(adata, use_highly_variable=False, random_state=RANDOM_SEED)
    adata.obsm["X_pca_original"] = adata.obsm["X_pca"].copy()
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)

    # ── 6. Leiden at all resolutions (+ optional precomputed reload) ─────
    resolutions = _resolution_range(RES_SCAN_MIN, RES_SCAN_MAX, RES_SCAN_STEP)
    leiden_keys = [f"leiden_{str(r).replace('.', '_')}" for r in resolutions]

    def _compute_all():
        for res in resolutions:
            key = f"leiden_{str(res).replace('.', '_')}"
            sc.tl.leiden(adata, resolution=res, key_added=key, random_state=RANDOM_SEED)
            log.info("  res=%.1f → %d clusters", res, adata.obs[key].nunique())

    if USE_PRECOMPUTED:
        ext_meta = (os.path.join(PRECOMPUTED_DIR, f"metadata_{sample_id}.tsv")
                    if PRECOMPUTED_DIR else "")
        if ext_meta and os.path.isfile(ext_meta):
            meta_source = ext_meta
            log.info("Reloading precomputed clusters from EXTERNAL: %s", meta_source)
        elif os.path.isfile(out_metadata):
            meta_source = out_metadata
            log.info("Reloading precomputed clusters from output: %s", meta_source)
        else:
            meta_source = None
            log.warning("use_precomputed is true but no metadata file found. Computing fresh.")

        if meta_source:
            saved_meta = pd.read_csv(meta_source, sep="\t", index_col=0, comment="#")
            for key in leiden_keys:
                if key in saved_meta.columns:
                    adata.obs[key] = pd.Categorical(
                        saved_meta[key].reindex(adata.obs_names).astype(int))
                    log.info("  Loaded %s: %d clusters", key, adata.obs[key].nunique())
                else:
                    log.warning("  %s not in saved metadata, computing …", key)
                    res = float(key.replace("leiden_", "").replace("_", "."))
                    sc.tl.leiden(adata, resolution=res, key_added=key, random_state=RANDOM_SEED)
        else:
            _compute_all()
    else:
        log.info("Computing Leiden for %d resolutions: %s", len(resolutions), list(resolutions))
        _compute_all()

    # ── 7. Save adata + metadata TSV ─────────────────────────────────────
    Path(out_adata).parent.mkdir(parents=True, exist_ok=True)
    log.info("Saving adata → %s", out_adata)
    adata.write(out_adata)

    meta_cols = ["sample", "region_annotation"] + [k for k in leiden_keys if k in adata.obs.columns]
    meta_df = adata.obs[meta_cols].copy()
    with open(out_metadata, "w") as f:
        f.write(f"# Generated: {datetime.now().isoformat()}\n")
        f.write(f"# Sample: {sample_id}\n")
        f.write(f"# Scanpy: {sc.__version__}\n")
        f.write(f"# Random seed: {RANDOM_SEED}\n")
        f.write(f"# n_cells: {adata.n_obs}\n")
        meta_df.to_csv(f, sep="\t")
    log.info("Preprocessing complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
