"""
generate_qupath_mer.py — MERSCOPE head: RGB morphology TIFF for QuPath + pixel<->micron
metadata, read straight from the Vizgen mosaic OME-TIFFs (no zarr).
=======================================================================================
Analogous to generate_qupath_x5k, but the source is the mosaic image folder rather than a
SpatialData zarr. DAPI + PolyT are always present; Cellbound1-3 are added when the Cell
Boundary Stain Kit was used. Emits three artefacts:

  * <sample>_morphology.tiff              RGB composite for QuPath annotation
  * <sample>_morphology_scalefactors.json pixel<->micron mapping (affine + scale + p0)
  * <sample>_background.npy               grey composite embedded VERBATIM as the contract
                                          background by prepare_input_mer

The grey background is written here (not rebuilt in prepare_input_mer) so the mosaic is
read/downsampled ONCE and the QuPath image, the contract background and the coordinate
metadata provably share a single grid. The coordinate math is likewise done once here and
only read back downstream — never re-derived — which is where MERSCOPE's affine offset is
kept from leaking into the region join.
"""
import json
import logging
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import tifffile

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
log = logging.getLogger("generate_qupath_mer")

# make the shared helpers importable (mirrors prepare_input_vhd's guard)
try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _here = os.getcwd()
sys.path.insert(0, _here)

try:
    from merscope_io import (find_merscope_files, read_micron_to_mosaic,
                             native_um_per_px, mosaic_to_micron_params, load_mosaic_stack)
    # generic (channel-agnostic) image helpers reused from the Xenium head
    from xenium_image import (normalize_channels, composite_rgb, rgb_uint8,
                              grey_projection_uint8)

    merscope_dir = str(snakemake.params.merscope_dir)
    out_tiff     = str(snakemake.output.qupath_image)
    out_meta     = str(snakemake.output.qupath_meta)
    out_bg       = str(snakemake.output.background)
    sample_id    = snakemake.params.sample_id
    z_index      = int(snakemake.params.z_index)
    target_um    = float(snakemake.params.hires_pixel_size_um)
    channels     = list(snakemake.params.channels)

    log.info("=" * 70)
    log.info("MERSCOPE QuPath composite for sample: %s", sample_id)
    log.info("=" * 70)

    files = find_merscope_files(merscope_dir)
    log.info("images dir: %s", files["images_dir"])

    # ── 1. Transform + downsample factor ─────────────────────────────────
    M = read_micron_to_mosaic(files["transform"])
    native = native_um_per_px(M)
    downsample = max(1, int(round(target_um / native)))
    log.info("Native mosaic %.4f µm/px; target %.3f µm/px → downsample x%d",
             native, target_um, downsample)

    # ── 2. Load + downsample the mosaic channels ─────────────────────────
    stack, names = load_mosaic_stack(files["images_dir"], z_index, channels, downsample)
    log.info("Loaded channels %s → stack %s", names, stack.shape)
    norm = normalize_channels(stack)

    # ── 3. RGB composite for QuPath (name-aware) + grey background ────────
    rgb_f, _mapping = composite_rgb(norm, names)
    rgb = rgb_uint8(rgb_f)
    Path(out_tiff).parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(out_tiff, rgb, compression="jpeg")
    log.info("Saved QuPath TIFF %s (%s)", out_tiff, rgb.shape)

    grey = grey_projection_uint8(norm)            # channel-agnostic; contract background
    np.save(out_bg, grey)
    log.info("Saved grey background %s (%s)", out_bg, grey.shape)

    # ── 4. Coordinate metadata (computed ONCE; read back by prepare_input) ─
    p0, pixel_to_um, affine6, scalef = mosaic_to_micron_params(M, downsample)
    meta = {
        "sample": sample_id,
        "z_index": z_index,
        "downsample": int(downsample),
        "qupath_image_height": int(rgb.shape[0]),
        "qupath_image_width": int(rgb.shape[1]),
        "channels": names,
        # image-origin frame (subtract p0 from raw centroids); see merscope_io
        "p0_micron": [float(p0[0]), float(p0[1])],
        "pixel_size_um": pixel_to_um,
        "qupath_pixel_to_um": pixel_to_um,             # x5k-compatible key (mean scale)
        "qupath_to_micron_affine": affine6,            # QuPath px -> image-origin micron
        "tissue_hires_scalef": scalef,                 # image-origin micron -> px
    }
    with open(out_meta, "w") as f:
        json.dump(meta, f, indent=2)
    log.info("Scale: 1 QuPath px = %.4f µm (p0 = %.2f, %.2f µm)",
             pixel_to_um, p0[0], p0[1])

except Exception:
    log.error("FAILED for %s:\n%s",
              getattr(snakemake.params, "sample_id", "?"), traceback.format_exc())
    raise
