"""
prepare_input_ate.py — Atera Headblock: build the standardized UNFILTERED CONTRACT
h5ad from the SpatialData zarr + QuPath GeoJSON.
==================================================================================
Analogous to prepare_input_x5k (Xenium 5K). Produces ONLY the contract:
  X                     raw integer counts (NO filtering / normalisation / clustering)
  obsm['spatial']       cell centroids (µm)
  uns['spatial'][sample] {images:{hires}, scalefactors}  (registered H&E when available,
                        otherwise the greyscale morphology composite)
  obs                   sample, cell_id, region_annotation

QC, normalisation, PCA/UMAP and clustering are the Coreblock's job (preprocess_umap) —
this script deliberately does none of them.

Region annotation accepts EITHER of the two QuPath images the head emits. Both live on
the morphology pixel grid, so polygons map to microns with a single scalar (no affine):

  {sample}_he_background.geojson  drawn on the H&E-on-morphology-grid background (the
                                  embedded contract background); px → µm = 1/hires_scalef
  {sample}_morphology.geojson     drawn on the morphology composite; px → µm from its
                                  scale-factor JSON

The H&E background is preferred when both exist, because it is the easier image to read.

The embedded contract background likewise prefers the registered H&E, which is far more
informative histologically than the fluorescence composite. It is used in the form
generate_qupath_he_ate resampled ONTO THE MORPHOLOGY GRID: sc.pl.spatial positions cells
with a single scalar factor and no rotation, and the two frames differ by ~90°, so the
raw H&E would misplace every cell. The morphology composite remains the fallback
whenever the H&E is not configured, is unreadable, or does not match the expected grid.
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
import pandas as pd

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
log = logging.getLogger("prepare_input_ate")

try:
    import spatialdata as sd

    zarr_path      = snakemake.params.zarr_path
    out_h5ad       = str(snakemake.output.h5ad)
    sample_id      = snakemake.params.sample_id
    geojson_dir    = snakemake.params.geojson_dir
    hires_level    = int(snakemake.params.hires_pyramid_level)
    pixel_size_um  = float(snakemake.params.pixel_size_um)
    qupath_meta    = getattr(snakemake.input, "qupath_meta", None)
    he_background  = getattr(snakemake.input, "he_background", None)

    log.info("=" * 70)
    log.info("Atera contract for sample: %s", sample_id)
    log.info("=" * 70)

    # ── 1. Read zarr ─────────────────────────────────────────────────────
    sdata = sd.read_zarr(zarr_path)
    log.info("SpatialData loaded from %s", zarr_path)

    # ── 2. Extract adata (RAW counts) ────────────────────────────────────
    adata = sdata.tables["table"].copy()
    adata.obs["sample"] = sample_id
    if "cell_id" not in adata.obs.columns:
        adata.obs["cell_id"] = adata.obs_names.astype(str)
    # cell_id must equal obs_names as strings — downstream joins key on it.
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    log.info("AnnData: %d cells × %d genes (raw, unfiltered)", adata.n_obs, adata.n_vars)

    # ── 2b. Contract dtype hygiene ───────────────────────────────────────
    # Newer pandas/anndata can emit extension dtypes (StringDtype index, nullable
    # Boolean/Int/Float) that break scanpy/scipy downstream (e.g. var_names.str
    # returning a BooleanArray inside calculate_qc_metrics). Normalise once via the
    # shared helper so the written contract is clean; the core readers apply the same
    # helper on load, so decoupled contracts and other heads are covered too.
    from contract_io import normalize_contract_dtypes
    adata = normalize_contract_dtypes(adata, log)

    # ── 3. Embed greyscale morphology composite in uns['spatial'] ────────
    #    Same channel composite as the QuPath TIFF, converted to greyscale, so the
    #    Coreblock's sc.pl.spatial has a background image. (Multi-channel data + masks
    #    stay in the zarr for spatialdata-based exploration outside the pipeline.)
    if "morphology_focus" in sdata.images:
        from spatialdata.models import get_channel_names
        img = sdata.images["morphology_focus"]
        ch_names = list(get_channel_names(img))
        scales = list(img.keys())
        level_idx = min(hires_level, len(scales) - 1)
        ds = img[scales[level_idx]]
        arr = ds[list(ds.data_vars)[0]].values
        if hasattr(arr, "compute"):
            arr = arr.compute()
        ds0 = img[scales[0]]
        full_shape = ds0[list(ds0.data_vars)[0]].shape           # (C, H, W)

        n_ch, h, w = arr.shape
        # Channel-AGNOSTIC greyscale morphology (mean projection of normalised channels);
        # robust to channel identity/order — see xenium_image.py. The multi-colour version
        # is only produced for the QuPath TIFF, not the data-bearing contract.
        from xenium_image import normalize_channels, grey_projection_uint8
        background = grey_projection_uint8(normalize_channels(arr))
        background_source = "morphology"
        # obsm['spatial'] is in MICRONS, so the scalefactor must convert micron -> level pixel:
        #   level_px = micron / pixel_size_um * (level_H / full_H)
        # (using only the pyramid factor would shrink the cells onto a corner of the image).
        hires_scalef = (h / full_shape[1]) / pixel_size_um

        # Prefer the registered H&E when the optional H&E rule ran: it is far more
        # informative histologically than the fluorescence composite. It has already been
        # RESAMPLED ONTO THIS MORPHOLOGY GRID by generate_qupath_he_ate — the raw H&E
        # cannot be used directly, because sc.pl.spatial has no rotation term and the two
        # frames differ by ~90°. The scalefactor is therefore unchanged.
        if he_background and os.path.isfile(str(he_background)):
            try:
                import tifffile
                he_img = tifffile.imread(str(he_background))
                if he_img.ndim == 3 and he_img.shape[:2] == (h, w):
                    background = he_img
                    background_source = "he"
                    log.info("Embedded registered H&E background (%s)", he_img.shape)
                else:
                    log.warning(
                        "H&E background shape %s does not match the morphology level "
                        "grid %s — keeping the morphology composite. Check that "
                        "atera.hires_pyramid_level is the same for both rules.",
                        getattr(he_img, "shape", None), (h, w))
                    del he_img
            except Exception as e:                                 # noqa: BLE001
                log.warning("Could not read the H&E background %s (%s) — keeping the "
                            "morphology composite.", he_background, e)
        elif he_background:
            log.warning("H&E background %s is missing or empty — keeping the morphology "
                        "composite.", he_background)

        adata.uns["spatial"] = {
            sample_id: {
                "images": {"hires": background},
                "scalefactors": {
                    "tissue_hires_scalef": float(hires_scalef),
                    "spot_diameter_fullres": 10.0,
                    "pixel_size_um": pixel_size_um,
                },
            }
        }
        adata.uns["morphology_channels"] = list(ch_names)
        adata.uns["background_image_source"] = background_source
        log.info("Embedded %s background (%s, scalef=%.6f)",
                 background_source, background.shape, hires_scalef)
        del arr, background
    else:
        log.warning("No morphology_focus in zarr — contract will have no embedded image "
                    "(the Coreblock will scatter on obsm['spatial']).")

    del sdata

    # ── 4. QuPath GeoJSON → region_annotation (optional) ─────────────────
    #    Prefer the H&E annotation when the optional H&E rule ran and the user drew on
    #    it; otherwise fall back to the morphology annotation.
    adata.obs["region_annotation"] = "Unlabeled"

    def _region_name(x):
        """QuPath classification -> region name, robust to how the GeoJSON reader hands
        back the nested object: a dict (geopandas/pyogrio) or a JSON string (older fiona)."""
        if isinstance(x, dict):
            return x.get("name", "Unlabeled")
        if isinstance(x, str):
            try:
                return json.loads(x).get("name", "Unlabeled")
            except (ValueError, AttributeError):
                return "Unlabeled"
        return "Unlabeled"

    bg_geojson = os.path.join(geojson_dir or "", f"{sample_id}_he_background.geojson")
    morph_geojson = os.path.join(geojson_dir or "", f"{sample_id}_morphology.geojson")

    # Micron<->pixel scale of the embedded background (set in section 3). Polygons drawn
    # on the contract-background image live on exactly this grid, so a single scalar maps
    # them to microns — no affine, which is why it is the most robust surface to annotate.
    _sf = adata.uns.get("spatial", {}).get(sample_id, {}).get("scalefactors", {})
    hires_scalef = _sf.get("tissue_hires_scalef")

    geojson_file, meta_file, space = None, None, None
    if geojson_dir and os.path.isfile(bg_geojson):
        # Drawn on {sample}_he_background.tiff — morphology-grid px at hires_pyramid_level.
        # Preferred: aligned with the cells by construction, needs only the scalefactor.
        geojson_file, meta_file, space = bg_geojson, None, "contract_background"
    elif geojson_dir and os.path.isfile(morph_geojson):
        geojson_file, meta_file, space = morph_geojson, str(qupath_meta or ""), "morphology"
    _present = [g for g in (bg_geojson, morph_geojson) if os.path.isfile(g)]
    if len(_present) > 1 and geojson_file:
        log.info("Multiple GeoJSONs present (%s) — using %s (space=%s).",
                 [os.path.basename(g) for g in _present],
                 os.path.basename(geojson_file), space)

    if geojson_file:
        log.info("Reading QuPath GeoJSON (%s space): %s", space, geojson_file)
        try:
            import geopandas as gpd
            from shapely.affinity import scale as shapely_scale

            meta = {}
            if space == "morphology":
                if meta_file and os.path.isfile(meta_file):
                    with open(meta_file) as f:
                        meta = json.load(f)
                else:
                    log.warning("Scale-factor JSON %s missing — falling back to 1 px = 1 µm, "
                                "which is almost certainly wrong.", meta_file)

            annot_gdf = gpd.read_file(geojson_file)
            if "classification" not in annot_gdf.columns:
                raise ValueError(
                    f"{geojson_file} has no 'classification' property — were the QuPath "
                    f"annotations assigned a class? Columns: {list(annot_gdf.columns)}")
            annot_gdf["region"] = annot_gdf["classification"].apply(_region_name)

            if space == "contract_background":
                # Polygons are on the embedded background's grid; the inverse of the
                # micron->pixel scalefactor takes them straight back to microns.
                if not hires_scalef:
                    raise ValueError(
                        "Annotated the contract background, but the contract has no "
                        "embedded scalefactor (no morphology_focus in the zarr?). Cannot "
                        "map polygons to microns.")
                px_to_um = 1.0 / float(hires_scalef)
                log.info("  Contract-background scale: 1 px = %.4f µm "
                         "(1 / tissue_hires_scalef)", px_to_um)
                annot_gdf["geometry"] = annot_gdf["geometry"].apply(
                    lambda g: shapely_scale(g, xfact=px_to_um, yfact=px_to_um, origin=(0, 0)))
            else:
                px_to_um = meta.get("qupath_pixel_to_um", 1.0)
                log.info("  Morphology scale: 1 QuPath px = %.4f µm", px_to_um)
                annot_gdf["geometry"] = annot_gdf["geometry"].apply(
                    lambda g: shapely_scale(g, xfact=px_to_um, yfact=px_to_um,
                                            origin=(0, 0)))

            spatial = adata.obsm["spatial"]
            cell_gdf = gpd.GeoDataFrame(
                {"cell_id": adata.obs["cell_id"].values},
                geometry=gpd.points_from_xy(spatial[:, 0], spatial[:, 1]))
            joined = gpd.sjoin(cell_gdf, annot_gdf[["region", "geometry"]],
                               how="left", predicate="within")
            joined = joined.drop_duplicates(subset="cell_id", keep="first")
            adata.obs["region_annotation"] = (
                joined.set_index("cell_id")["region"]
                .reindex(adata.obs["cell_id"].values, fill_value="Unlabeled").values)
            counts = adata.obs["region_annotation"].value_counts().to_dict()
            log.info("Regions: %s", counts)
            # A transform error usually shows up as "everything is Unlabeled" rather
            # than as an exception, so say so explicitly.
            if set(counts) == {"Unlabeled"}:
                log.warning("No cell fell inside any annotated polygon. If regions were "
                            "drawn, the %s coordinate transform is suspect — check the "
                            "scale factors in %s.", space, meta_file)
            adata.uns["region_annotation_source"] = {
                "geojson": geojson_file, "space": space,
            }
        except Exception as e:
            log.warning("QuPath annotation failed (%s); region_annotation = Unlabeled", e)
            adata.obs["region_annotation"] = "Unlabeled"
    else:
        log.info("No GeoJSON at %s or %s — region_annotation = Unlabeled "
                 "(regions optional)", bg_geojson, morph_geojson)

    # ── 5. Write the unfiltered contract (X = raw counts) ────────────────
    Path(out_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)
    log.info("Wrote contract %s  (%d cells × %d genes, raw counts)",
             out_h5ad, adata.n_obs, adata.n_vars)

except Exception:
    log.error("FAILED for %s:\n%s", snakemake.params.sample_id, traceback.format_exc())
    raise
