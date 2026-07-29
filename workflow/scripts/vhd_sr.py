"""Space Ranger output discovery, shared by the two Visium HD head steps
(generate_qupath_vhd + prepare_input_vhd) so both discover the SR hires image the
same way — one definition, no drift.
"""
import os
from glob import glob

# hires tissue image: the fixed Space Ranger locations, then a recursive catch-all
# for unusual layouts (first match wins). Both head steps use this exact list.
HIRES_IMAGE_PATTERNS = (
    "outs/segmented_outputs/spatial/tissue_hires_image.png",
    "outs/spatial/tissue_hires_image.png",
    "segmented_outputs/spatial/tissue_hires_image.png",
    "spatial/tissue_hires_image.png",
    "**/tissue_hires_image.png",
)


def find_sr_file(sr_outdir, *patterns):
    """Return the first file under sr_outdir matching any glob pattern (tried in order).

    recursive=True so a ``**`` catch-all pattern works; it is a no-op for the plain
    patterns, so all callers are safe.
    """
    for pattern in patterns:
        matches = glob(os.path.join(sr_outdir, pattern), recursive=True)
        if matches:
            return matches[0]
    tried = ", ".join(patterns)
    raise FileNotFoundError(f"Could not find any of [{tried}] under {sr_outdir}")
