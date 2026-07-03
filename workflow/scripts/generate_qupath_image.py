"""
generate_qupath_image.py — Visium HD HEAD (step 2).
===================================================
Copies the Space Ranger hires tissue image out of the (version-varying) SR tree
to a stable path, so it can be opened and annotated in QuPath. The exported
geojson (named "{sample}_tissue_hires_image.geojson") is later consumed by
prepare_input. Pure file copy — no data processing.
"""
import logging
import os
import shutil
import sys
import traceback
from glob import glob
from pathlib import Path

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
log = logging.getLogger("generate_qupath_image")


def find_sr_file(sr_outdir, *patterns):
    for pattern in patterns:
        matches = glob(os.path.join(sr_outdir, pattern))
        if matches:
            return matches[0]
    raise FileNotFoundError(f"None of {patterns} found under {sr_outdir}")


try:
    sample_id = snakemake.params.sample_id
    sr_outdir = snakemake.params.sr_outdir
    out_image = str(snakemake.output.qupath_image)

    src = find_sr_file(
        sr_outdir,
        "outs/segmented_outputs/spatial/tissue_hires_image.png",
        "outs/spatial/tissue_hires_image.png",
        "segmented_outputs/spatial/tissue_hires_image.png",
        "spatial/tissue_hires_image.png",
    )
    Path(out_image).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, out_image)
    log.info("QuPath image for %s: %s -> %s", sample_id, src, out_image)

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
