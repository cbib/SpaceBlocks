"""
validate_input.py — Input-contract validation for the modular CORE.
===================================================================
Runs once per sample, BEFORE qc_sweep / run_upstream. It asserts that the
unfiltered, segmented h5ad produced by the (technology-specific) HEAD satisfies
the core contract, and records whether an image is embedded so that spatial
plots downstream can choose between an image-backed `sc.pl.spatial` and a plain
coordinate scatter. It NEVER modifies the data — file preparation is a head
concern; this rule only checks.

Outcome
-------
* Hard requirements violated  → raises (no report written → Snakemake aborts the
  whole DAG before any analysis runs).
* All hard requirements met   → writes a JSON report with the assessed metadata
  (most importantly ``image_present`` / ``image_mode``) that downstream rules can
  read; soft issues are logged and recorded as ``warnings``.

Contract checked
----------------
HARD : loads as AnnData; >0 cells and >0 genes; X present, finite, non-negative
       and (optionally) integer-valued raw counts; obsm[spatial_key] present,
       (n_obs, >=2), finite; sample_key in obs; obs_names & var_names unique;
       if uns[spatial_key] exists it must be a well-formed image (images +
       scalefactors per library), else it is a contract violation.
SOFT : region_annotation absent (regions are optional); no mito genes (pct_mt
       unavailable); no embedded image (→ scatter fallback recorded, not an error).
"""
import json
import logging
import sys
import traceback
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp

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
log = logging.getLogger("validate_input")


def _x_min(X):
    """Minimum of X, handling sparse and dense."""
    if sp.issparse(X):
        return float(X.min()) if X.nnz else 0.0
    return float(np.asarray(X).min()) if np.size(X) else 0.0


def _is_integer_counts(X, max_check=200_000):
    """Are the (non-zero) values integer-valued? Sampled for speed."""
    d = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if d.size == 0:
        return True
    if not np.all(np.isfinite(d)):
        return False
    if d.size > max_check:
        idx = np.random.default_rng(0).choice(d.size, max_check, replace=False)
        d = d[idx]
    return bool(np.allclose(d, np.round(d)))


def _valid_image_entry(entry):
    """A library entry is a usable image iff it carries images + scalefactors."""
    return (isinstance(entry, dict)
            and isinstance(entry.get("images"), dict) and len(entry["images"]) > 0
            and isinstance(entry.get("scalefactors"), dict)
            and len(entry["scalefactors"]) > 0)


try:
    h5ad_path   = str(snakemake.input.h5ad)
    sample_id   = str(snakemake.params.sample_id)
    sample_key  = str(snakemake.params.sample_key)
    spatial_key = str(snakemake.params.spatial_key)
    require_region = bool(snakemake.params.require_region)
    require_raw    = bool(snakemake.params.require_raw_counts)
    mito_prefix = tuple(snakemake.params.mito_prefix)
    out_report  = str(snakemake.output.report)

    log.info("=" * 70)
    log.info("Validating input contract: sample=%s", sample_id)
    log.info("  file=%s", h5ad_path)
    log.info("=" * 70)

    errors, warnings = [], []

    # ── Load ─────────────────────────────────────────────────────────────
    adata = sc.read_h5ad(h5ad_path)
    n_obs, n_vars = int(adata.n_obs), int(adata.n_vars)
    if n_obs == 0:
        errors.append("0 cells in object")
    if n_vars == 0:
        errors.append("0 genes in object")

    # ── X = raw counts ───────────────────────────────────────────────────
    if adata.X is None:
        errors.append("X is None (expected raw counts)")
    else:
        if _x_min(adata.X) < 0:
            errors.append("X has negative values (expected non-negative raw counts)")
        if require_raw and not _is_integer_counts(adata.X):
            errors.append("X is not integer-valued — looks normalized, not raw "
                          "counts (set contract.require_raw_counts=false to allow)")

    # ── Spatial coordinates (segmented cells in space) ───────────────────
    has_coords = spatial_key in adata.obsm
    if not has_coords:
        errors.append(f"obsm['{spatial_key}'] missing (need spatial coordinates)")
    else:
        coords = np.asarray(adata.obsm[spatial_key])
        if coords.ndim != 2 or coords.shape[1] < 2:
            errors.append(f"obsm['{spatial_key}'] is not (n_obs, >=2)")
        elif coords.shape[0] != n_obs:
            errors.append(f"obsm['{spatial_key}'] row count != n_obs")
        elif not np.isfinite(coords).all():
            errors.append(f"obsm['{spatial_key}'] contains non-finite coordinates")

    # ── Required obs / uniqueness ────────────────────────────────────────
    if sample_key not in adata.obs.columns:
        errors.append(f"obs['{sample_key}'] missing (sample identifier required)")
    if not adata.obs_names.is_unique:
        errors.append("obs_names are not unique")
    if not adata.var_names.is_unique:
        errors.append("var_names are not unique")

    # ── Region annotation (optional) ─────────────────────────────────────
    region_present = "region_annotation" in adata.obs.columns
    if require_region and not region_present:
        errors.append("obs['region_annotation'] required but missing")
    elif not region_present:
        warnings.append("obs['region_annotation'] absent — regions optional, "
                        "region-based plots will be skipped downstream")

    # ── Mito genes (optional metric) ─────────────────────────────────────
    mito_available = bool(adata.var_names.str.startswith(mito_prefix).any())
    if not mito_available:
        warnings.append(f"no genes match mito prefix {mito_prefix} — "
                        "pct_counts_mt will be unavailable")

    # ── Image presence + structure (assessed, not prepared) ──────────────
    image_present, image_mode, image_library_ids = False, "scatter", []
    uns_sp = adata.uns.get(spatial_key, None)
    if isinstance(uns_sp, dict) and len(uns_sp) > 0:
        libs = list(uns_sp.keys())
        if all(_valid_image_entry(uns_sp[lib]) for lib in libs):
            image_present, image_mode, image_library_ids = True, "embedded", libs
        else:
            errors.append(f"uns['{spatial_key}'] present but malformed: every "
                          "library needs non-empty 'images' and 'scalefactors' "
                          "(fix in the head, or omit the image for scatter mode)")
    else:
        warnings.append("no embedded image — spatial plots will scatter on "
                        f"obsm['{spatial_key}'] (scatter mode)")

    # ── Verdict ──────────────────────────────────────────────────────────
    for w in warnings:
        log.warning("  %s", w)
    passed = len(errors) == 0

    report = {
        "sample": sample_id,
        "passed": passed,
        "file": h5ad_path,
        "n_obs": n_obs,
        "n_vars": n_vars,
        "spatial_coords": has_coords,
        "spatial_key": spatial_key,
        "image_present": image_present,
        "image_mode": image_mode,             # "embedded" → sc.pl.spatial w/ image
        "image_library_ids": image_library_ids,  # "scatter"  → scatter on coords
        "sample_key": sample_key,
        "region_annotation_present": region_present,
        "mito_available": mito_available,
        "errors": errors,
        "warnings": warnings,
    }

    if not passed:
        for e in errors:
            log.error("  CONTRACT VIOLATION: %s", e)
        # Do NOT write the gate output — Snakemake then fails the DAG here.
        raise ValueError(
            f"Input contract validation FAILED for {sample_id}: "
            + "; ".join(errors))

    Path(out_report).parent.mkdir(parents=True, exist_ok=True)
    with open(out_report, "w") as fh:
        json.dump(report, fh, indent=2)
    log.info("Validation PASSED for %s  (image_present=%s, mode=%s)",
             sample_id, image_present, image_mode)

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
