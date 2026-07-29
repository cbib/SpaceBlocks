"""
convert_zarr_x5k.py — Convert Xenium output bundle to SpatialData zarr.
=======================================================================
Uses spatialdata_io.xenium() exactly as in the 10x 5k analysis notebook.
The zarr store becomes the canonical data object containing:
  - morphology_focus (MultiscaleSpatialImage, all channels)
  - cell_labels, nucleus_labels (segmentation masks)
  - cell_boundaries, nucleus_boundaries (polygon shapes)
  - table (AnnData with count matrix + spatial coords)
"""
import logging, os, sys
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
log = logging.getLogger("convert_zarr")

import spatialdata_io as sio

xenium_dir = str(snakemake.input.xenium_dir)
zarr_path = snakemake.params.zarr_path
sample_id = snakemake.params.sample_id

log.info("Converting Xenium → zarr for sample %s", sample_id)
log.info("  Input:  %s", xenium_dir)
log.info("  Output: %s", zarr_path)

sdata = sio.xenium(
    xenium_dir,
    cells_boundaries=True,
    nucleus_boundaries=True,
    cells_as_circles=True,
    cells_labels=True,
    nucleus_labels=True,
    transcripts=False,          # large; not needed for pipeline
    morphology_mip=False,       # we use morphology_focus only
    morphology_focus=True,
    aligned_images=True,
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
