"""Validate an optional region-annotation GeoJSON without modifying it."""
import json
import logging
import sys
import traceback
from pathlib import Path

import geopandas as gpd

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
log = logging.getLogger("process_geojson")

# ── Parameters ───────────────────────────────────────────────────────────────
in_geojson = str(snakemake.input.geojson)
sample_id = str(snakemake.params.sample_id)
out_report = str(snakemake.output.report)

try:
    log.info("=" * 70)
    log.info("Processing region GeoJSON for sample: %s", sample_id)
    log.info("  input: %s", in_geojson)
    log.info("=" * 70)

    # ── 1. Load ────────────────────────────────────────────────────────
    gdf = gpd.read_file(in_geojson)
    n_total = len(gdf)
    log.info("Loaded %d feature(s)", n_total)

    if n_total == 0:
        raise ValueError(f"'{in_geojson}' contains no features.")

    # ── 2. Validate ─────────────────────────────────────────────────────
    if "classification" not in gdf.columns:
        raise ValueError(
            f"'{in_geojson}' has no 'classification' property. Assign a QuPath "
            f"class to the annotations before exporting. Columns: {list(gdf.columns)}"
        )

    missing_mask = gdf.geometry.isna()
    n_missing = int(missing_mask.sum())
    empty_mask = gdf.geometry.is_empty & ~missing_mask
    n_empty = int(empty_mask.sum())
    invalid_mask = ~gdf.geometry.is_valid & ~missing_mask & ~empty_mask
    n_invalid = int(invalid_mask.sum())
    if n_missing or n_empty or n_invalid:
        raise ValueError(
            f"'{in_geojson}' contains unusable geometries: missing={n_missing}, "
            f"empty={n_empty}, invalid={n_invalid}. Fix them in QuPath and re-export; "
            "the workflow does not rewrite source annotations."
        )

    # ── 3. Write success report ─────────────────────────────────────────
    report = {
        "sample": sample_id,
        "passed": True,
        "file": str(Path(in_geojson).resolve()),
        "feature_count": n_total,
        "geometry_types": sorted(gdf.geometry.geom_type.unique().tolist()),
        "classification_present": True,
    }
    Path(out_report).parent.mkdir(parents=True, exist_ok=True)
    Path(out_report).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    log.info("GeoJSON validation passed; report written to %s", out_report)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
