"""
xenium_image.py — robust morphology-channel compositing for the Xenium5k head.
==============================================================================
Shared by export_qupath_x5k (RGB TIFF for QuPath annotation) and prepare_input_x5k
(greyscale background embedded in the contract), so the two can never diverge.

This module:

  * normalises each channel independently (1–99th percentile),
  * assigns a colour by channel NAME (nuclear→blue, membrane/boundary→red,
    RNA/interior→green), with a positional fallback palette when names are
    unrecognised, handling ANY channel count,
  * logs the channel→role→colour mapping so it can be verified from the run log,
  * and for the CONTRACT background uses a channel-AGNOSTIC mean projection, so the
    data-bearing artefact never depends on the colour assumption at all.

None of this affects counts, coordinates, or region_annotation — the composite is
only a visualisation / annotation aid — but it makes the image safe across setups.
"""
import logging

import numpy as np

log = logging.getLogger("xenium_image")

# Colour roles → RGB. Nuclear = blue, membrane/boundary = red-magenta, RNA = green,
# cytoplasm/interior protein = amber.
_ROLE_RGB = {
    "nuclear":   (0.0, 0.0, 1.0),
    "membrane":  (1.0, 0.0, 0.4),
    "rna":       (0.0, 0.9, 0.3),
    "cytoplasm": (0.95, 0.7, 0.0),
    "other":     (0.6, 0.6, 0.6),
}
# Distinct fallback colours when channel names are not recognised at all.
_FALLBACK = [(0, 0, 1), (0, 0.9, 0.3), (1, 0, 0.4), (0.9, 0.9, 0),
             (0, 0.8, 0.9), (0.9, 0, 0.9)]


def _role_for(name):
    """Map a channel name to a colour role (case-insensitive substring match), or None.
    Covers the standard Xenium multimodal-segmentation stack, e.g.
    ['DAPI', 'ATP1A1/CD45/E-Cadherin', '18S', 'AlphaSMA/Vimentin']."""
    n = str(name).lower()
    if any(k in n for k in ("dapi", "nucl", "hoechst")):
        return "nuclear"
    if any(k in n for k in ("18s", "rna", "interior")):
        return "rna"
    if any(k in n for k in ("boundary", "membrane", "atp1a1", "cadherin",
                            "e-cad", "ecad", "cd45")):
        return "membrane"
    if any(k in n for k in ("alphasma", "sma", "vimentin", "actin", "cytoplasm")):
        return "cytoplasm"
    return None


def normalize_channels(arr):
    """Percentile-normalise each channel of a (C, H, W) stack to [0, 1] float32."""
    arr = np.asarray(arr)
    out = np.zeros(arr.shape, dtype=np.float32)
    for i in range(arr.shape[0]):
        ch = arr[i].astype(np.float32)
        if ch.max() > 0:
            p1, p99 = np.percentile(ch[ch > 0], [1, 99])
            out[i] = np.clip((ch - p1) / max(p99 - p1, 1), 0, 1)
    return out


def composite_rgb(norm, channel_names=None):
    """Name-aware RGB composite from a normalised (C, H, W) stack.

    Returns (rgb float32 HxWx3 in [0, 1], mapping list of (name, role, rgb)). Logs the
    mapping. Unrecognised names fall back to a distinct-colour palette (if NO channel
    is recognised) or grey (if some are), so nothing is silently dropped."""
    n_ch, h, w = norm.shape
    names = list(channel_names) if channel_names is not None else []
    names += [f"ch{i}" for i in range(len(names), n_ch)]

    named_any = any(_role_for(names[i]) is not None for i in range(n_ch))
    rgb = np.zeros((h, w, 3), dtype=np.float32)
    mapping = []
    for i in range(n_ch):
        role = _role_for(names[i])
        if role is not None:
            col = _ROLE_RGB[role]
        elif not named_any:
            col = _FALLBACK[i % len(_FALLBACK)]
            role = "positional"
        else:
            col = _ROLE_RGB["other"]
            role = "other"
        for c in range(3):
            rgb[:, :, c] += norm[i] * col[c]
        mapping.append((names[i], role, col))

    log.info("Morphology composite: %d channel(s), name-aware=%s", n_ch, named_any)
    for nm, role, col in mapping:
        log.info("  '%s' → %-10s rgb=%s", nm, role, tuple(round(x, 2) for x in col))
    return np.clip(rgb, 0, 1), mapping


def rgb_uint8(rgb):
    """(H,W,3) float[0,1] → uint8 for tifffile."""
    return (np.clip(rgb, 0, 1) * 255).astype(np.uint8)


def grey_projection_uint8(norm):
    """Channel-AGNOSTIC mean projection of a normalised (C, H, W) stack → (H, W, 3)
    uint8. Used for the contract background so it never depends on channel identity."""
    grey = norm.mean(axis=0) if norm.shape[0] > 0 else np.zeros(norm.shape[1:])
    return (np.stack([grey, grey, grey], axis=-1) * 255).astype(np.uint8)
