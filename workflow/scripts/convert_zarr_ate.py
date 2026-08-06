"""
convert_zarr_ate.py — Convert an Atera output bundle to SpatialData zarr.
=========================================================================
Atera is built on the Xenium platform and its `outs/` bundle follows the Xenium
Onboard Analysis v4 layout, so `spatialdata_io.xenium()` reads it directly. The zarr
store becomes the canonical data object containing:
  - morphology_focus (MultiscaleSpatialImage, all channels, named from the OME XML)
  - cell_labels, nucleus_labels (segmentation masks)
  - cell_boundaries, nucleus_boundaries (polygon shapes)
  - table (AnnData with count matrix + spatial coords)

Two Atera-specific notes:

  * `aligned_images=False`. The registered H&E ships as a SEPARATE download rather
    than inside the bundle, so there is nothing for the reader to auto-discover. It
    is handled by generate_qupath_he_ate, which keeps the multi-GB image out of the
    zarr entirely.

  * The morphology channel files are named chNNNN_<stain>.ome.tif rather than
    morphology_focus_NNNN.ome.tif. That is the XOA v4 layout and is handled natively
    by spatialdata-io >= 0.7; the preflight below fails with an actionable message
    rather than a reader traceback when an older version is in the environment.
"""
import json, logging, sys
from pathlib import Path
import shutil

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
                    handlers=log_handlers)
log = logging.getLogger("convert_zarr_ate")

import spatialdata_io as sio

atera_dir = str(snakemake.input.atera_dir)
zarr_path = snakemake.params.zarr_path
sample_id = snakemake.params.sample_id

log.info("Converting Atera → zarr for sample %s", sample_id)
log.info("  Input:  %s", atera_dir)
log.info("  Output: %s", zarr_path)

# ── Preflight: report the bundle's provenance before the expensive read ──────
# `analysis_sw_version` decides which format branch the reader takes; logging it
# makes format regressions obvious in the run log rather than at the traceback.
specs_path = Path(atera_dir) / "experiment.xenium"
if specs_path.is_file():
    try:
        specs = json.loads(specs_path.read_text())
        log.info("  Bundle: chemistry=%s  panel=%s (%s targets)  analysis_sw=%s",
                 specs.get("chemistry_version"), specs.get("panel_name"),
                 specs.get("panel_num_targets_predesigned"),
                 specs.get("analysis_sw_version"))
        log.info("  Reported pixel_size=%s µm, %s cells",
                 specs.get("pixel_size"), specs.get("num_cells"))
    except Exception as e:                                     # noqa: BLE001
        log.warning("Could not parse %s (%s) — continuing.", specs_path, e)
else:
    log.warning("No experiment.xenium at %s — is atera_dir pointing at the outs/ "
                "directory of the bundle?", atera_dir)

# The XOA v4 morphology layout needs spatialdata-io >= 0.7. Fail here with an
# actionable message instead of deep inside the reader.
focus_dir = Path(atera_dir) / "morphology_focus"
if focus_dir.is_dir():
    focus_files = sorted(p.name for p in focus_dir.glob("*.ome.tif"))
    log.info("  morphology_focus/: %s", focus_files)
    if focus_files and not any(f.startswith("morphology_focus_") for f in focus_files):
        try:
            from importlib.metadata import version as _pkg_version
            import packaging.version as _pv
            sio_version = _pkg_version("spatialdata-io")
            if _pv.parse(sio_version) < _pv.parse("0.7"):
                sys.exit(
                    f"[convert_zarr_ate] This bundle uses the XOA v4 morphology layout "
                    f"({focus_files[0]}), which requires spatialdata-io >= 0.7; the "
                    f"environment has {sio_version}. Update workflow/envs/atera.yaml.")
            log.info("  XOA v4 morphology layout, spatialdata-io %s — supported.",
                     sio_version)
        except Exception as e:                                 # noqa: BLE001
            log.warning("Could not verify the spatialdata-io version (%s).", e)

sdata = sio.xenium(
    atera_dir,
    cells_boundaries=True,
    nucleus_boundaries=True,
    cells_as_circles=True,
    cells_labels=True,
    nucleus_labels=True,
    transcripts=False,          # ~10 GB parquet; not needed for cell-level analysis
    morphology_mip=False,       # we use morphology_focus only
    morphology_focus=True,
    aligned_images=False,       # the H&E is a separate download — see module docstring
    cells_table=True,
)
log.info("SpatialData loaded:\n%s", sdata)

# Delete morphology_mip if present (saves disk)
if "morphology_mip" in sdata.images:
    del sdata.images["morphology_mip"]
    log.info("Removed morphology_mip to save disk space")

zarr_path_obj = Path(zarr_path)
zarr_path_obj.parent.mkdir(parents=True, exist_ok=True)

log.info("Path exists before write? %s", zarr_path_obj.exists())

if zarr_path_obj.exists():
    log.warning("Removing existing path: %s", zarr_path)
    shutil.rmtree(zarr_path_obj)

sdata.write(zarr_path)

log.info("Zarr written to %s", zarr_path)

# Completion marker — Snakemake tracks THIS (touched only after a successful write),
# so a job that dies mid-write leaves no marker and is correctly re-run.
Path(snakemake.output.done).touch()
log.info("Wrote completion marker %s", snakemake.output.done)
