"""
generate_qupath_vhd.py — Visium HD HEAD (step 2)
==================================================
Copies the Space Ranger hires tissue image out of the (version-varying) SR output
tree to a stable, predictable path so it can be opened and annotated in QuPath.
The manually exported "{sample}_tissue_hires_image.geojson" is later consumed by
prepare_input to build the region_annotation contract column.
"""

import logging
import os
import shutil
import sys
import traceback
from glob import glob
from pathlib import Path

# ── Logging (mirror the other head/core scripts) ─────────────────────────────
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


# Space Ranger file discovery, shared with the other VHD head step (no drift).
try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:                      # very old Snakemake
    _here = os.getcwd()
sys.path.insert(0, _here)
from vhd_sr import find_sr_file, HIRES_IMAGE_PATTERNS


try:
    sample_id = str(snakemake.params.sample_id)
    sr_outdir = str(snakemake.params.sr_outdir)
    out_png   = str(snakemake.output.qupath_image)

    log.info("=" * 70)
    log.info("QuPath image copy: sample=%s", sample_id)
    log.info("  SR outdir: %s", sr_outdir)
    log.info("=" * 70)

    # Same fallbacks as prepare_input, plus a recursive catch-all for unusual
    # SR layouts (first match wins).
    hires_src = find_sr_file(sr_outdir, *HIRES_IMAGE_PATTERNS)
    log.info("  found hires image: %s", hires_src)

    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(hires_src, out_png)
    log.info("  copied → %s", out_png)
    log.info("QuPath image ready for %s", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s",
              getattr(snakemake.params, "sample_id", "?"), traceback.format_exc())
    raise
