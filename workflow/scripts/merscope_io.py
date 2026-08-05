"""
merscope_io.py — shared MERSCOPE (Vizgen MERFISH) file discovery + coordinate helpers.
=======================================================================================
Imported by BOTH MERSCOPE head steps — generate_qupath_mer (QuPath composite + the
pixel<->micron metadata) and prepare_input_mer (contract build) — so the two can never
drift on which files they read or how mosaic pixels map to microns. Same role vhd_sr.py
plays for the Visium HD head.

MERSCOPE region layout (one region == one sample):

    <region>/
      <prefix>cell_by_gene<suffix>.csv        transcripts-per-cell matrix (raw counts)
      <prefix>cell_metadata<suffix>.csv       per-cell centroids (center_x/y, µm) + QC
      images/
        micron_to_mosaic_pixel_transform.csv  3x3 affine, micron -> mosaic pixel
        mosaic_DAPI_z{0..6}.tif               nuclear stain   (ALWAYS present)
        mosaic_PolyT_z{0..6}.tif              total-RNA stain (ALWAYS present)
        mosaic_Cellbound{1,2,3}_z*.tif        membrane stains (only with the boundary kit)

Coordinate convention used by the head
--------------------------------------
The mosaic transform is an affine micron -> mosaic-pixel with a NON-zero translation
(mosaic pixel (0,0) is not micron (0,0)). To keep the QuPath round-trip identical to the
Xenium head — a bare isotropic scale with origin (0,0), no offset — the head expresses
cell coordinates in *microns measured from mosaic pixel (0,0)* (the "image-origin frame",
obtained by subtracting `p0` below). In that frame the micron->pixel map loses its
translation, so the embedded background needs only a scalefactor and the GeoJSON round-trip
is a pure linear map. The transform math is done ONCE here / in generate_qupath_mer and
written to the QuPath meta JSON; prepare_input_mer only reads it back.
"""
import glob
import logging
import os
import re

import numpy as np
import pandas as pd

log = logging.getLogger("merscope_io")

# DAPI + PolyT are always imaged; Cellbound1-3 appear only when the Cell Boundary Stain
# Kit was used. This is also the default channel order for the QuPath RGB composite.
DEFAULT_CHANNELS = ["DAPI", "PolyT", "Cellbound1", "Cellbound2", "Cellbound3"]


def _one(root, *patterns, required=True, what=""):
    """Return the single file under `root` matching the first pattern that hits."""
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(root, pat)))
        if hits:
            if len(hits) > 1:
                log.warning("Multiple matches for %s (%s); using %s", what, pat, hits[0])
            return hits[0]
    if required:
        raise FileNotFoundError(
            f"No {what or 'file'} under {root} matching any of: {patterns}")
    return None


def find_merscope_files(root):
    """Resolve the MERSCOPE inputs under a region directory.

    Filenames carry a dataset-specific prefix/suffix in public releases (e.g.
    'datasets_..._cell_by_gene_S1R1.csv'), so we glob rather than hard-code names.
    Returns a dict with keys: counts, meta, transform, images_dir.
    """
    images_dir = _one(root, "images", "*/images", required=True, what="images dir")
    return {
        "counts": _one(root, "*cell_by_gene*.csv", "*/*cell_by_gene*.csv",
                       what="cell_by_gene matrix"),
        "meta": _one(root, "*cell_metadata*.csv", "*/*cell_metadata*.csv",
                     what="cell_metadata"),
        "transform": _one(images_dir, "*micron_to_mosaic_pixel_transform.csv",
                          what="micron_to_mosaic transform"),
        "images_dir": images_dir,
    }


def read_micron_to_mosaic(transform_csv):
    """Read the 3x3 micron->mosaic-pixel affine M.

    [px_x, px_y, 1]^T = M @ [µm_x, µm_y, 1]^T  (px_x = image column, px_y = image row).
    Whitespace- or comma-delimited; three rows of three values.
    """
    try:
        M = np.loadtxt(transform_csv)
    except ValueError:
        M = np.loadtxt(transform_csv, delimiter=",")
    return np.asarray(M, dtype=float).reshape(3, 3)


def native_um_per_px(M):
    """Full-resolution mosaic pixel size in microns, from the affine's linear part."""
    return float(1.0 / np.sqrt(abs(np.linalg.det(M[:2, :2]))))


def mosaic_to_micron_params(M, downsample):
    """Derive the head's coordinate mapping from the affine M and an integer image
    `downsample` factor.

    Returns:
      p0          (2,) micron coordinate of mosaic pixel (0,0). Subtract from raw
                  centroids to move them into the image-origin frame.
      pixel_to_um float, output-grid µm per pixel (area-preserving mean).
      affine6     [a, b, d, e, xoff, yoff] mapping a downsampled-image pixel (x=col,
                  y=row) -> image-origin micron, for geopandas .affine_transform().
                  xoff = yoff = 0 by construction (no offset in the image-origin frame),
                  and for the usual diagonal transform this reduces to the Xenium head's
                  `geometry.scale(pixel_to_um, pixel_to_um, origin=(0, 0))`.
      scalef      float, image-origin micron -> downsampled pixel (tissue_hires_scalef).
    """
    lin = M[:2, :2]              # micron -> pixel linear part
    t = M[:2, 2]                 # translation
    inv = np.linalg.inv(lin)     # pixel -> micron linear part
    p0 = -inv @ t                # micron at mosaic pixel (0,0)
    f = float(downsample)
    # downsampled pixel q -> full pixel (q*f) -> image-origin micron:  µm' = f * inv @ q
    a, b = f * inv[0, 0], f * inv[0, 1]
    d, e = f * inv[1, 0], f * inv[1, 1]
    affine6 = [float(a), float(b), float(d), float(e), 0.0, 0.0]
    pixel_to_um = float(np.sqrt(abs(a * e - b * d)))
    scalef = float(1.0 / pixel_to_um)
    return p0.astype(float), pixel_to_um, affine6, scalef


def downsample_mean(arr, f):
    """Block-mean downsample a 2-D array by integer factor f (dependency-free)."""
    f = int(f)
    if f <= 1:
        return np.asarray(arr, dtype=np.float32)
    h = arr.shape[0] - arr.shape[0] % f
    w = arr.shape[1] - arr.shape[1] % f
    a = np.asarray(arr[:h, :w], dtype=np.float32)
    return a.reshape(h // f, f, w // f, f).mean(axis=(1, 3))


def load_mosaic_stack(images_dir, z_index, channels, downsample):
    """Load requested mosaic channels at z=z_index, block-mean downsampled by `downsample`.

    Channels absent on disk (e.g. Cellbound without the boundary kit) are skipped, so a
    run stained with DAPI+PolyT only still yields a background. Planes are read one at a
    time and the raw plane freed immediately, to bound peak memory on whole-slide mosaics.
    Returns (stack (C,H,W) float32, names).
    """
    import tifffile
    planes, names = [], []
    for ch in channels:
        # tolerate a dataset prefix on public releases (e.g. 'datasets_..._mosaic_DAPI_z3.tif')
        hits = sorted(glob.glob(os.path.join(images_dir, f"*mosaic_{ch}_z{z_index}.tif")))
        if not hits:
            log.info("  channel %s absent at z%d — skipping", ch, z_index)
            continue
        raw = tifffile.imread(hits[0])
        planes.append(downsample_mean(raw, downsample))
        names.append(ch)
        del raw
    if not planes:
        raise FileNotFoundError(
            f"No mosaic channels {channels} at z{z_index} under {images_dir}")
    # after crop-to-multiple, planes can differ by a pixel; trim to the common shape
    h = min(p.shape[0] for p in planes)
    w = min(p.shape[1] for p in planes)
    stack = np.stack([p[:h, :w] for p in planes], axis=0).astype(np.float32)
    return stack, names


def read_cell_table(counts_csv, meta_csv, drop_blank=True):
    """Read the Vizgen cell-by-gene matrix + cell metadata into aligned frames.

    Returns (counts DataFrame cells×genes, meta DataFrame cells×fields), indices cast to
    str and intersected. `Blank-*` control-probe columns are dropped from counts by default
    (they measure the false-positive rate and are not real genes).
    """
    counts = pd.read_csv(counts_csv, index_col=0)
    meta = pd.read_csv(meta_csv, index_col=0)
    counts.index = counts.index.astype(str)
    meta.index = meta.index.astype(str)
    if drop_blank:
        blanks = [c for c in counts.columns if str(c).lower().startswith("blank")]
        if blanks:
            log.info("Dropping %d Blank control columns", len(blanks))
            counts = counts.drop(columns=blanks)
    common = counts.index.intersection(meta.index)
    if len(common) == 0:
        raise RuntimeError("No overlapping cell IDs between counts and metadata")
    return counts.loc[common], meta.loc[common]
