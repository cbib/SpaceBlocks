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
out_report     = str(snakemake.output.report)
SAMPLE_META    = dict(getattr(snakemake.params, "sample_meta", {}) or {})
output_dir     = str(Path(out_adata).parent)

# ── Optional per-sample QC threshold overrides ───────────────────────────────
# Start from the analysis.* config defaults (already read above); override any of
# them for THIS sample from the optional thresholds TSV (first column = sample
# name; columns = min_counts/min_genes/min_cells/max_counts/max_pct_mt). Missing
# samples or blank cells fall back to the config default, key-by-key.
_thresholds = {
    "min_counts": MIN_COUNTS,
    "min_cells":  MIN_CELLS,
    "min_genes":  MIN_GENES,
    "max_counts": snakemake.params.max_counts,   # config default (may be None)
    "max_pct_mt": snakemake.params.max_pct_mt,   # config default (may be None)
}
THRESHOLDS_SOURCE = "analysis.* config defaults"
_th_tsv = getattr(snakemake.input, "thresholds_tsv", None)
if _th_tsv:
    try:
        _tdf = pd.read_csv(str(_th_tsv), sep="\t", comment="#")
        _scol = _tdf.columns[0]
        _tdf[_scol] = _tdf[_scol].astype(str)
        _row = _tdf[_tdf[_scol] == str(sample_id)]
        if len(_row):
            _row = _row.iloc[0]
            _applied = []
            for _k in ("min_counts", "min_cells", "min_genes", "max_counts", "max_pct_mt"):
                if _k in _tdf.columns and pd.notna(_row[_k]) and str(_row[_k]).strip() != "":
                    _thresholds[_k] = _row[_k]
                    _applied.append(_k)
            THRESHOLDS_SOURCE = (f"{_th_tsv} (row for {sample_id}: "
                                 f"{', '.join(_applied) if _applied else 'no columns'})")
            log.info("Per-sample threshold overrides for %s: %s",
                     sample_id, _applied or "none")
        else:
            log.info("No override row for %s in %s; using config defaults",
                     sample_id, _th_tsv)
    except Exception as e:
        log.warning("Could not read per-sample thresholds (%s): %s; using config defaults",
                    _th_tsv, e)


def _opt_int(v):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else int(float(v))


def _opt_float(v):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)


# Effective, resolved cut-offs used for filtering (min_* required, max_* optional)
MIN_COUNTS = _opt_int(_thresholds["min_counts"])
MIN_CELLS  = _opt_int(_thresholds["min_cells"])
MIN_GENES  = _opt_int(_thresholds["min_genes"])
MAX_COUNTS = _opt_int(_thresholds["max_counts"])
MAX_PCT_MT = _opt_float(_thresholds["max_pct_mt"])


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

    # Stamp the core sample-sheet metadata columns into obs (one constant value
    # per sample). These propagate through annotation, integration, and pseudobulk
    # so downstream rules can group / colour by them (see config.extra_annotations).
    for _dk, _dv in SAMPLE_META.items():
        adata.obs[_dk] = pd.Categorical([str(_dv)] * adata.n_obs)
    if SAMPLE_META:
        log.info("  stamped sample-sheet columns into obs: %s", list(SAMPLE_META))

    # ── 2. QC metrics + filter ───────────────────────────────────────────
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hb"], inplace=True, log1p=False)

    # Pre-filter snapshot for the per-sample report (observed feature ranges are
    # taken on the full, unfiltered dataset).
    mito_present   = bool(adata.var["mt"].any())
    n_cells_before = int(adata.n_obs)
    n_genes_before = int(adata.n_vars)

    def _range(series):
        v = np.asarray(series.values, dtype=float)
        return (float(np.min(v)), float(np.max(v))) if v.size else (float("nan"), float("nan"))

    rng_counts = _range(adata.obs["total_counts"])
    rng_genes  = _range(adata.obs["n_genes_by_counts"])
    rng_mt     = (_range(adata.obs["pct_counts_mt"])
                  if "pct_counts_mt" in adata.obs.columns else (float("nan"), float("nan")))
    rng_ncells = (_range(adata.var["n_cells_by_counts"])
                  if "n_cells_by_counts" in adata.var.columns else (float("nan"), float("nan")))

    log.info("QC (source: %s): min_counts=%s min_cells=%s min_genes=%s "
             "max_counts=%s max_pct_mt=%s", THRESHOLDS_SOURCE,
             MIN_COUNTS, MIN_CELLS, MIN_GENES, MAX_COUNTS, MAX_PCT_MT)
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    if MAX_COUNTS is not None:
        sc.pp.filter_cells(adata, max_counts=MAX_COUNTS)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    if MAX_PCT_MT is not None:
        if mito_present:
            _n_pre = int(adata.n_obs)
            adata = adata[adata.obs["pct_counts_mt"] <= MAX_PCT_MT].copy()
            log.info("  max_pct_mt=%.3f removed %d cells", MAX_PCT_MT, _n_pre - adata.n_obs)
        else:
            log.warning("  max_pct_mt set but no MT- genes present; mito filter skipped")
    log.info("After QC: %d cells, %d genes", adata.n_obs, adata.n_vars)

    n_cells_after = int(adata.n_obs)
    n_genes_after = int(adata.n_vars)

    # ── 3. Region annotation plot (already joined upstream) ──────────────
    if "region_annotation" not in adata.obs.columns:
        adata.obs["region_annotation"] = "Unlabeled"
    if isinstance(REGION_COLORS, dict) and REGION_COLORS:
        adata.obs["region_annotation"] = adata.obs["region_annotation"].astype("category")
        cats = adata.obs["region_annotation"].cat.categories
        adata.uns["region_annotation_colors"] = [REGION_COLORS.get(str(c), "#cccccc") for c in cats]
    try:
        # Always pass an explicit spot_size (cell centroids are far smaller than
        # the scalefactor-derived Visium spot size, so without this the dots render
        # invisibly and only the tissue image shows); add library_id to draw them
        # over the embedded image when one is present.
        spatial_kw = {"spot_size": 20}
        if isinstance(adata.uns.get("spatial"), dict) and adata.uns["spatial"]:
            spatial_kw["library_id"] = list(adata.uns["spatial"].keys())[0]
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

    # ── 8. Per-sample QC report ──────────────────────────────────────────
    # Cells before/after filtering, observed min/max of each QC feature (on the
    # unfiltered data), and the cut-offs actually applied (min_* required,
    # max_* only when set). Long format (sample, metric, value) so per-sample
    # reports concatenate trivially into a cohort table.
    def _fmt(x):
        try:
            return "NA" if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 3)
        except (TypeError, ValueError):
            return "NA"

    report_rows = [
        ("n_cells_before",         n_cells_before),
        ("n_cells_after",          n_cells_after),
        ("n_cells_removed",        n_cells_before - n_cells_after),
        ("n_genes_before",         n_genes_before),
        ("n_genes_after",          n_genes_after),
        ("total_counts_min",       _fmt(rng_counts[0])),
        ("total_counts_max",       _fmt(rng_counts[1])),
        ("n_genes_by_counts_min",  _fmt(rng_genes[0])),
        ("n_genes_by_counts_max",  _fmt(rng_genes[1])),
        ("pct_counts_mt_min",      _fmt(rng_mt[0])),
        ("pct_counts_mt_max",      _fmt(rng_mt[1])),
        ("n_cells_per_gene_min",   _fmt(rng_ncells[0])),
        ("n_cells_per_gene_max",   _fmt(rng_ncells[1])),
        ("cutoff_min_counts",      MIN_COUNTS if MIN_COUNTS is not None else "NA"),
        ("cutoff_max_counts",      MAX_COUNTS if MAX_COUNTS is not None else "NA"),
        ("cutoff_min_genes",       MIN_GENES  if MIN_GENES  is not None else "NA"),
        ("cutoff_min_cells",       MIN_CELLS  if MIN_CELLS  is not None else "NA"),
        ("cutoff_max_pct_mt",      MAX_PCT_MT if MAX_PCT_MT is not None else "NA"),
        ("thresholds_source",      THRESHOLDS_SOURCE),
    ]
    # Record the sample-sheet metadata columns parsed for this sample.
    report_rows += [(f"sheet:{_dk}", _dv) for _dk, _dv in SAMPLE_META.items()]
    pd.DataFrame([{"sample": sample_id, "metric": m, "value": v}
                  for m, v in report_rows]).to_csv(out_report, sep="\t", index=False)
    log.info("Wrote per-sample report → %s", out_report)

    log.info("Preprocessing complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
