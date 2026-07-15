"""
export_qupath_image.py — Extract composite TIFF from zarr for QuPath.

Reads morphology channels at a configurable pyramid level, composites
into an RGB TIFF, and saves a scale-factor JSON for coordinate mapping
between QuPath pixel coords and the micron coordinate system.
"""
import json, logging, os, sys
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
log = logging.getLogger("export_qupath")

import spatialdata as sd
from spatialdata.models import get_channel_names

zarr_path = snakemake.params.zarr_path
out_tiff = str(snakemake.output.qupath_image)
out_meta = str(snakemake.output.qupath_meta)
qupath_level = int(snakemake.params.qupath_pyramid_level)
sample_id = snakemake.params.sample_id

log.info("Exporting QuPath image for %s at pyramid level %d", sample_id, qupath_level)

sdata = sd.read_zarr(zarr_path)

if "morphology_focus" not in sdata.images:
    log.error("No morphology_focus in zarr — creating empty placeholder")
    Path(out_tiff).parent.mkdir(parents=True, exist_ok=True)
    Path(out_tiff).touch()
    with open(out_meta, "w") as f:
        json.dump({"error": "no morphology_focus found"}, f)
    sys.exit(0)

img = sdata.images["morphology_focus"]
ch_names = list(get_channel_names(img))
scales = list(img.keys())
pixel_size_um = 0.2125  # Xenium default

# Read at requested pyramid level
level_idx = min(qupath_level, len(scales) - 1)
ds = img[scales[level_idx]]
var_name = list(ds.data_vars)[0]
arr = ds[var_name].values
if hasattr(arr, "compute"):
    arr = arr.compute()
log.info("  Level %d: shape %s, channels: %s", level_idx, arr.shape, ch_names)

# Full-res shape for scale factor
ds0 = img[scales[0]]
full_shape = ds0[list(ds0.data_vars)[0]].shape  # (C, H, W)

# Composite RGB
n_ch, h, w = arr.shape
rgb = np.zeros((h, w, 3), dtype=np.float32)

for i in range(n_ch):
    ch = arr[i].astype(np.float32)
    if ch.max() > 0:
        p1, p99 = np.percentile(ch[ch > 0], [1, 99])
        ch = np.clip((ch - p1) / max(p99 - p1, 1), 0, 1)
    else:
        ch[:] = 0
    if i == 0:      # DAPI → blue
        rgb[:, :, 2] += ch
    elif i == 1:    # Interior protein → green
        rgb[:, :, 1] += ch * 0.8
    elif i == 2:    # Boundary → magenta
        rgb[:, :, 0] += ch * 0.95
        rgb[:, :, 2] += ch * 0.4
    elif i == 3:    # Interior RNA → yellow
        rgb[:, :, 0] += ch * 0.6
        rgb[:, :, 1] += ch * 0.6

rgb = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)

Path(out_tiff).parent.mkdir(parents=True, exist_ok=True)
tifffile.imwrite(out_tiff, rgb, compression="jpeg")
log.info("  Saved %s (%s)", out_tiff, rgb.shape)

# Scale metadata
qupath_to_fullres = full_shape[1] / h
meta = {
    "sample": sample_id,
    "pyramid_level": level_idx,
    "qupath_image_height": h,
    "qupath_image_width": w,
    "fullres_height": int(full_shape[1]),
    "fullres_width": int(full_shape[2]),
    "qupath_to_fullres_scale": float(qupath_to_fullres),
    "pixel_size_um": pixel_size_um,
    "qupath_pixel_to_um": float(qupath_to_fullres * pixel_size_um),
    "channels": ch_names,
}
with open(out_meta, "w") as f:
    json.dump(meta, f, indent=2)
log.info("  Scale: 1 QuPath px = %.4f µm", meta["qupath_pixel_to_um"])
