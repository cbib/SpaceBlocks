"""
process_geojson.py – CORE: load, validate, and stage a sample's region GeoJSON.
================================================================================
Reads the OPTIONAL region-annotation GeoJSON for a sample (either
{sample}_tissue_hires_image.geojson or {sample}_morphology.geojson, resolved
upstream by _find_geojson in common.smk — this script only ever sees a single,
already-resolved path, since presence/naming is validated before the rule is even
scheduled). Sanity-checks the geometries, tags each feature with the sample ID,
and writes it into the sample's output directory for downstream ROI-based
analyses to consume.
"""
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
in_geojson  = str(snakemake.input.geojson)
sample_id   = snakemake.params.sample_id
out_geojson = str(snakemake.output.geojson)

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

    # ── 2. Validate / repair geometries ───────────────────────────────────
    invalid_mask = ~gdf.geometry.is_valid
    n_invalid = int(invalid_mask.sum())
    if n_invalid:
        log.warning("%d/%d invalid geometries found — attempting buffer(0) repair",
                    n_invalid, n_total)
        gdf.loc[invalid_mask, "geometry"] = gdf.loc[invalid_mask, "geometry"].buffer(0)
        still_invalid = int((~gdf.geometry.is_valid).sum())
        if still_invalid:
            raise ValueError(
                f"{still_invalid}/{n_total} geometries remain invalid after repair "
                f"— inspect '{in_geojson}'.")

    # ── 3. Tag with sample of origin ─────────────────────────────────────
    # Survives later concatenation across samples without losing provenance.
    gdf["sample_id"] = sample_id

    # ── 4. Save ────────────────────────────────────────────────────────
    Path(out_geojson).parent.mkdir(parents=True, exist_ok=True)
    log.info("Saving %d feature(s) → %s", n_total, out_geojson)
    gdf.to_file(out_geojson, driver="GeoJSON")

    log.info("GeoJSON processing complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise