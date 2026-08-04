"""
atera_he.py — H&E ↔ morphology registration helpers for the Atera HeadBlock.
============================================================================
Shared by generate_qupath_he_ate (writes the composed transform + keypoint QA into
the scale-factor JSON) and prepare_input_ate (applies it to QuPath polygons), so the
two can never diverge — the same relationship xenium_image.py has for the morphology
composite.

The 10x `*_he_alignment.csv` is a row-major 3x3 homogeneous matrix mapping
FULL-RESOLUTION H&E pixels -> FULL-RESOLUTION morphology (DAPI) pixels:

    [ m00 m01 m02 ]   [x_he]     [x_morph]
    [ m10 m11 m12 ] . [y_he]  =  [y_morph]
    [  0   0   1  ]   [  1 ]     [   1   ]

For Atera the 2x2 block is a similarity transform (scale + rotation, no shear).
Direction was verified against `*_keypoints.csv`, whose `alignmentX/Y` columns are
H&E pixels and `fixedX/Y` columns are morphology pixels — `keypoint_residuals()`
below re-runs that check at export time so a future change of convention fails
loudly instead of silently displacing every region annotation.

None of this touches counts or coordinates; it only positions the region polygons.
"""
import logging

import numpy as np

log = logging.getLogger("atera_he")


def read_alignment_matrix(path):
    """Read a 10x `*_he_alignment.csv` (3x3, no header) into a (3, 3) float array."""
    M = np.loadtxt(str(path), delimiter=",", dtype=float)
    if M.shape != (3, 3):
        raise ValueError(
            f"Expected a 3x3 alignment matrix in {path}, got shape {M.shape}. "
            "This should be the 10x '<sample>_he_alignment.csv' supplemental file.")
    if not np.allclose(M[2], [0.0, 0.0, 1.0], atol=1e-6):
        log.warning("Alignment matrix bottom row is %s, expected [0, 0, 1] — the file "
                    "may not be a plain affine; continuing with the 2x3 block.", M[2])
    return M


def similarity_params(M):
    """(scale, rotation_degrees) of the 2x2 block, for logging and sanity checks.

    Meaningful only when the block is a similarity transform (Atera's is). The
    rotation sign follows image coordinates (y increasing downwards)."""
    a, b = float(M[0, 0]), float(M[0, 1])
    scale = float(np.hypot(a, b))
    rotation = float(np.degrees(np.arctan2(b, a)))
    return scale, rotation


def is_similarity(M, rtol=1e-3):
    """True when the 2x2 block is [[a, b], [-b, a]] (scale + rotation, no shear)."""
    return (np.isclose(M[0, 0], M[1, 1], rtol=rtol)
            and np.isclose(M[0, 1], -M[1, 0], rtol=rtol))


def he_to_morphology(M, xy):
    """Apply the alignment matrix to an (N, 2) array of full-res H&E pixel coords.
    Returns an (N, 2) array of full-res morphology pixel coords."""
    xy = np.asarray(xy, dtype=float).reshape(-1, 2)
    ones = np.ones((xy.shape[0], 1))
    return (np.hstack([xy, ones]) @ M.T)[:, :2]


def warp_he_to_morphology_grid(he_rgb, M, he_level_scale, morph_level_scale, target_hw):
    """Resample a downsampled H&E RGB image onto the morphology pixel grid.

    `sc.pl.spatial` places cells over the background with a single scalar factor and no
    rotation, so an H&E used as a contract background must first be brought into the
    morphology frame — otherwise the ~90° rotation between the two puts every cell in
    the wrong place.

    Parameters
    ----------
    he_rgb            (H, W, 3) uint8 H&E read at some pyramid level.
    M                 (3, 3) alignment matrix: full-res H&E px -> full-res morphology px.
    he_level_scale    he_rgb height / full-resolution H&E height.
    morph_level_scale target grid height / full-resolution morphology height.
    target_hw         (H, W) of the morphology pyramid level to resample onto.

    Returns an (H, W, 3) uint8 array on the target grid, white outside the H&E.

    The composition, in (x, y): a target pixel maps back to full-res morphology by
    dividing by `morph_level_scale`, then through M^-1 to full-res H&E, then times
    `he_level_scale` to the loaded level. scipy's ``affine_transform`` samples the INPUT
    at ``matrix @ output_coord + offset`` in (row, col) order, so the (x, y) form is
    transposed by S = [[0, 1], [1, 0]] before being handed over.
    """
    from scipy.ndimage import affine_transform as _ndi_affine

    Minv = np.linalg.inv(np.asarray(M, dtype=float))
    A, t = Minv[:2, :2], Minv[:2, 2]

    B = (float(he_level_scale) / float(morph_level_scale)) * A   # (x, y)
    c = float(he_level_scale) * t                                # (x, y)

    S = np.array([[0.0, 1.0], [1.0, 0.0]])
    B_rc = S @ B @ S
    c_rc = np.array([c[1], c[0]])

    out = np.empty((int(target_hw[0]), int(target_hw[1]), 3), dtype=np.uint8)
    for ch in range(3):
        out[..., ch] = _ndi_affine(
            he_rgb[..., ch], B_rc, offset=c_rc, output_shape=tuple(int(v) for v in target_hw),
            order=1, mode="constant", cval=255,          # white = empty slide
        ).astype(np.uint8)
    log.info("Warped H&E %s -> morphology grid %s", he_rgb.shape[:2], tuple(target_hw))
    return out


def keypoint_residuals(M, keypoints_path, pixel_size_um):
    """Validate the alignment matrix against the 10x `*_keypoints.csv`.

    Pushes each (alignmentX, alignmentY) H&E point through M and compares with the
    paired (fixedX, fixedY) morphology point. Returns a dict of residual statistics
    in both pixels and microns, or None when the file is absent or unreadable.

    Hand-placed keypoints on a whole-slide scan normally land within a few pixels;
    anything at the scale of a tissue region means the convention has changed."""
    import pandas as pd

    try:
        kp = pd.read_csv(str(keypoints_path))
    except Exception as e:                                     # noqa: BLE001
        log.warning("Could not read keypoints %s (%s) — skipping alignment QA.",
                    keypoints_path, e)
        return None

    needed = {"fixedX", "fixedY", "alignmentX", "alignmentY"}
    if not needed.issubset(kp.columns):
        log.warning("Keypoints file %s lacks columns %s — skipping alignment QA.",
                    keypoints_path, sorted(needed - set(kp.columns)))
        return None
    if kp.empty:
        log.warning("Keypoints file %s has no rows — skipping alignment QA.",
                    keypoints_path)
        return None

    predicted = he_to_morphology(M, kp[["alignmentX", "alignmentY"]].to_numpy())
    observed = kp[["fixedX", "fixedY"]].to_numpy(dtype=float)
    dist_px = np.linalg.norm(predicted - observed, axis=1)

    return {
        "n_keypoints": int(len(kp)),
        "residual_px_mean": float(dist_px.mean()),
        "residual_px_max": float(dist_px.max()),
        "residual_um_mean": float(dist_px.mean() * pixel_size_um),
        "residual_um_max": float(dist_px.max() * pixel_size_um),
    }
