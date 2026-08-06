"""
generate_qupath_he_ate.py — Registered H&E → QuPath TIFF, contract background, transform.
=========================================================================================
The Atera supplemental files include a registered H&E whole-slide image and a 3x3 affine
(`*_he_alignment.csv`) mapping full-resolution H&E pixels onto full-resolution morphology
pixels. This rule:

  1. reads the H&E at a configurable pyramid level (the full-resolution OME-TIFF is tens
     of GB — it is never loaded whole, and never enters the zarr),
  2. writes an RGB TIFF the user can annotate in QuPath,
  3. composes QuPath-px → full-res H&E px → morphology px → µm into a single affine and
     stores it in the scale-factor JSON for prepare_input_ate to apply to polygons,
  4. resamples the H&E onto the MORPHOLOGY pixel grid, producing the background image
     prepare_input_ate embeds in the contract,
  5. validates the alignment against `*_keypoints.csv` and warns on large residuals.

Why step 4 exists
-----------------
`sc.pl.spatial` positions cells over the background with a single scalar scale factor —
there is no rotation term. The Atera H&E is related to the morphology frame by a
similarity transform that is very nearly a 90° rotation, so embedding the H&E as-is
would place every cell in the wrong location. Warping it onto the morphology grid makes
the H&E drop-in compatible with the same `tissue_hires_scalef` the morphology composite
uses, and keeps the contract's coordinate convention identical across HeadBlocks.

Why step 5 exists
-----------------
The alignment convention is not documented by 10x and was established empirically. If a
future release flips it, silent misplacement of every region annotation is the failure
mode. The keypoint check turns that into a loud warning at export time, before anyone
spends an afternoon in QuPath.

NOTE for users: the H&E and morphology QuPath TIFFs for one sample are related by roughly
a 90° rotation and will not look alike. That is expected — annotate whichever you prefer.
"""
import json, logging, sys
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
log = logging.getLogger("generate_qupath_he_ate")

from atera_he import (read_alignment_matrix, similarity_params, is_similarity,
                      keypoint_residuals, warp_he_to_morphology_grid)

he_image_path = str(snakemake.input.he_image)
he_align_path = str(snakemake.input.he_alignment)
he_keypoints = snakemake.params.he_keypoints
zarr_path = snakemake.params.zarr_path
out_meta = str(snakemake.output.qupath_meta)
out_background = str(snakemake.output.he_background)
he_level = int(snakemake.params.he_pyramid_level)
hires_level = int(snakemake.params.hires_pyramid_level)
pixel_size_um = float(snakemake.params.pixel_size_um)
residual_warn_px = float(snakemake.params.residual_warn_px)
sample_id = snakemake.params.sample_id

log.info("=" * 70)
log.info("Atera H&E QuPath image for sample: %s", sample_id)
log.info("=" * 70)

# ── 1. Alignment matrix ──────────────────────────────────────────────────────
M = read_alignment_matrix(he_align_path)
scale, rotation = similarity_params(M)
log.info("Alignment matrix from %s", he_align_path)
log.info("  H&E px → morphology px: scale=%.5f, rotation=%.2f°, similarity=%s",
         scale, rotation, is_similarity(M))
log.info("  Implied H&E resolution: %.4f µm/px", scale * pixel_size_um)
if not is_similarity(M):
    log.warning("The 2x2 block is not a similarity transform (shear present). The "
                "composed affine still applies, but check the registration.")

# ── 2. Read the H&E at a pyramid level ───────────────────────────────────────
# Never read level 0: the Atera H&E is a multi-GB whole-slide pyramid.
with tifffile.TiffFile(he_image_path) as tf:
    series = tf.series[0]
    levels = getattr(series, "levels", [series])
    n_levels = len(levels)
    level_idx = min(he_level, n_levels - 1)
    if he_level > level_idx:
        log.warning("Requested H&E pyramid level %d but only %d level(s) exist — "
                    "using level %d.", he_level, n_levels, level_idx)
    full_shape = tuple(levels[0].shape)
    log.info("H&E pyramid: %d level(s), full-res shape %s", n_levels, full_shape)
    arr = levels[level_idx].asarray()

arr = np.squeeze(arr)
# Normalise to (H, W, 3): whole-slide RGB may arrive channel-first depending on writer.
if arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[-1] not in (3, 4):
    arr = np.moveaxis(arr, 0, -1)
if arr.ndim == 3 and arr.shape[-1] == 4:
    arr = arr[..., :3]                      # drop alpha
if arr.ndim == 2:
    arr = np.stack([arr] * 3, axis=-1)      # greyscale scan → RGB
if arr.dtype != np.uint8:
    finite = arr[np.isfinite(arr)] if np.issubdtype(arr.dtype, np.floating) else arr
    hi = float(np.percentile(finite, 99.5)) if finite.size else 1.0
    arr = (np.clip(arr.astype(np.float32) / max(hi, 1e-6), 0, 1) * 255).astype(np.uint8)

h, w = arr.shape[0], arr.shape[1]
log.info("  Level %d: %d x %d px", level_idx, h, w)

Path(out_background).parent.mkdir(parents=True, exist_ok=True)

# ── 3. Compose QuPath px → µm ────────────────────────────────────────────────
# full_shape is the level-0 H&E shape; take its height however the axes are ordered.
_channel_first = full_shape[0] in (3, 4)
fullres_h = int(full_shape[1] if _channel_first else full_shape[0])
fullres_w = int(full_shape[2] if _channel_first else full_shape[1])
qupath_to_fullres = fullres_h / h

# ── 4. Contract background: resample the H&E onto the morphology grid ────────
# prepare_input_ate embeds this instead of the greyscale morphology composite. It must
# land on exactly the grid that rule computes its scalefactor from, so both read the
# morphology pyramid at `hires_pyramid_level`.
import spatialdata as sd

sdata = sd.read_zarr(zarr_path)
if "morphology_focus" not in sdata.images:
    log.error("No morphology_focus in the zarr — cannot define the target grid for the "
              "H&E background. Writing an empty placeholder; prepare_input_ate will fall "
              "back to scattering on obsm['spatial'].")
    Path(out_background).touch()
    m_level, target_hw = None, None
else:
    mimg = sdata.images["morphology_focus"]
    mscales = list(mimg.keys())
    m_level = min(hires_level, len(mscales) - 1)
    mds = mimg[mscales[m_level]]
    m_shape = mds[list(mds.data_vars)[0]].shape                 # (C, H, W)
    mds0 = mimg[mscales[0]]
    m_full = mds0[list(mds0.data_vars)[0]].shape                # (C, H, W)
    target_hw = (int(m_shape[1]), int(m_shape[2]))
    morph_level_scale = target_hw[0] / int(m_full[1])
del sdata

if target_hw is not None:
    he_level_scale = h / fullres_h
    log.info("Warping H&E onto the morphology grid: %s → %s (morphology level %d)",
             (h, w), target_hw, m_level)
    if he_level_scale < morph_level_scale:
        log.warning("The H&E level (%.4f of full-res) is coarser than the morphology "
                    "target level (%.4f) — the background will be upsampled and look "
                    "soft. Lower atera.he_pyramid_level for a sharper background.",
                    he_level_scale, morph_level_scale)
    warped = warp_he_to_morphology_grid(arr, M, he_level_scale, morph_level_scale,
                                        target_hw)
    tifffile.imwrite(out_background, warped, compression="jpeg")
    log.info("  Saved contract background %s (%s)", out_background, warped.shape)
    del warped
del arr

# ── 5. Keypoint QA ───────────────────────────────────────────────────────────
qa = None
if he_keypoints and Path(str(he_keypoints)).is_file():
    qa = keypoint_residuals(M, he_keypoints, pixel_size_um)
    if qa:
        log.info("Alignment QA over %d keypoints: mean %.1f px (%.2f µm), "
                 "max %.1f px (%.2f µm)",
                 qa["n_keypoints"], qa["residual_px_mean"], qa["residual_um_mean"],
                 qa["residual_px_max"], qa["residual_um_max"])
        if qa["residual_px_max"] > residual_warn_px:
            log.warning(
                "Max keypoint residual %.1f px exceeds the %.0f px threshold. Either "
                "the registration is poor for this sample, or the alignment-matrix "
                "convention has changed in a newer Atera release. Region annotations "
                "drawn on the H&E may be displaced, and the embedded background may be "
                "misregistered — verify before trusting them, or annotate the "
                "morphology image instead.",
                qa["residual_px_max"], residual_warn_px)
    else:
        log.warning("Keypoint QA could not be computed — the composed transform is "
                    "unverified for this sample.")
else:
    log.warning("No keypoints file configured or found (atera.he_keypoints) — the "
                "composed H&E transform will be applied UNVERIFIED. Supplying the "
                "10x '<sample>_keypoints.csv' is strongly recommended.")

# ── 6. Metadata ──────────────────────────────────────────────────────────────
meta = {
    "sample": sample_id,
    # prepare_input_ate branches on this: polygons drawn here need the affine below,
    # not the plain scalar scaling used for the morphology image.
    "annotation_space": "contract_background",
    "pyramid_level": level_idx,
    "qupath_image_height": int(h),
    "qupath_image_width": int(w),
    "fullres_height": fullres_h,
    "fullres_width": fullres_w,
    "qupath_to_fullres_scale": float(qupath_to_fullres),
    "pixel_size_um": pixel_size_um,
    "he_um_per_px_fullres": float(scale * pixel_size_um),
    "alignment_scale": float(scale),
    "alignment_rotation_deg": float(rotation),
    "alignment_matrix": [[float(v) for v in row] for row in M],
    # Ready to hand straight to shapely.affinity.affine_transform.
    # The morphology pyramid level the background was warped onto. prepare_input_ate
    # checks this against its own hires level before embedding.
    "background_morphology_level": (int(m_level) if m_level is not None else None),
    "background_height": (int(target_hw[0]) if target_hw else None),
    "background_width": (int(target_hw[1]) if target_hw else None),
    "keypoint_qa": qa,
}
with open(out_meta, "w") as f:
    json.dump(meta, f, indent=2)
log.info("Wrote %s", out_meta)
