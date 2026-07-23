"""
prepare_input_x5k.py — Xenium5k Headblock: build the standardized UNFILTERED CONTRACT
h5ad from the SpatialData zarr + QuPath GeoJSON.

Analogous to prepare_input.py (Visium HD). Produces ONLY the contract:
  X                     raw integer counts (NO filtering / normalisation / clustering)
  obsm['spatial']       cell centroids (µm)
  uns['spatial'][sample] {images:{hires}, scalefactors}  (greyscale morphology composite)
  obs                   sample, cell_id, region_annotation

QC, normalisation, PCA/UMAP and clustering are the Coreblock's job (preprocess_umap),
exactly as for Visium HD — this script deliberately does none of them.
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
log = logging.getLogger("prepare_input_x5k")

try:
    import spatialdata as sd

    zarr_path      = snakemake.params.zarr_path
    out_h5ad       = str(snakemake.output.h5ad)
    sample_id      = snakemake.params.sample_id
    geojson_dir    = snakemake.params.geojson_dir
    hires_level    = int(snakemake.params.hires_pyramid_level)
    pixel_size_um  = float(snakemake.params.pixel_size_um)
    qupath_meta    = getattr(snakemake.input, "qupath_meta", None)

    log.info("=" * 70)
    log.info("Xenium5k contract for sample: %s", sample_id)
    log.info("=" * 70)

    # ── 1. Read zarr ─────────────────────────────────────────────────────
    sdata = sd.read_zarr(zarr_path)
    log.info("SpatialData loaded from %s", zarr_path)

    # ── 2. Extract adata (RAW counts) ────────────────────────────────────
    adata = sdata.tables["table"].copy()
    adata.obs["sample"] = sample_id
    if "cell_id" not in adata.obs.columns:
        adata.obs["cell_id"] = adata.obs_names.astype(str)
    log.info("AnnData: %d cells × %d genes (raw, unfiltered)", adata.n_obs, adata.n_vars)

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
        grey_rgb = grey_projection_uint8(normalize_channels(arr))
        # obsm['spatial'] is in MICRONS, so the scalefactor must convert micron -> level pixel:
        #   level_px = micron / pixel_size_um * (level_H / full_H)
        # (using only the pyramid factor would shrink the cells onto a corner of the image).
        hires_scalef = (h / full_shape[1]) / pixel_size_um

        adata.uns["spatial"] = {
            sample_id: {
                "images": {"hires": grey_rgb},
                "scalefactors": {
                    "tissue_hires_scalef": float(hires_scalef),
                    "spot_diameter_fullres": 10.0,
                    "pixel_size_um": pixel_size_um,
                },
            }
        }
        adata.uns["morphology_channels"] = list(ch_names)
        log.info("Embedded greyscale morphology image (%s, scalef=%.6f)",
                 grey_rgb.shape, hires_scalef)
        del arr, grey_rgb
    else:
        log.warning("No morphology_focus in zarr — contract will have no embedded image "
                    "(the Coreblock will scatter on obsm['spatial']).")

    del sdata

    # ── 4. QuPath GeoJSON → region_annotation (optional) ─────────────────
    adata.obs["region_annotation"] = "Unlabeled"
    geojson_file = os.path.join(geojson_dir or "", f"{sample_id}_morphology.geojson")
    if geojson_dir and os.path.isfile(geojson_file):
        log.info("Reading QuPath GeoJSON: %s", geojson_file)
        try:
            import geopandas as gpd
            px_to_um = 1.0
            if qupath_meta and os.path.isfile(str(qupath_meta)):
                with open(str(qupath_meta)) as f:
                    px_to_um = json.load(f).get("qupath_pixel_to_um", 1.0)
            annot_gdf = gpd.read_file(geojson_file)
            annot_gdf["region"] = annot_gdf["classification"].apply(
                lambda x: json.loads(x).get("name", "Unlabeled")
                if isinstance(x, str) else "Unlabeled")
            annot_gdf["geometry"] = annot_gdf["geometry"].scale(
                xfact=px_to_um, yfact=px_to_um, origin=(0, 0))
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
            log.info("Regions: %s",
                     adata.obs["region_annotation"].value_counts().to_dict())
        except Exception as e:
            log.warning("QuPath annotation failed (%s); region_annotation = Unlabeled", e)
            adata.obs["region_annotation"] = "Unlabeled"
    else:
        log.info("No GeoJSON at %s — region_annotation = Unlabeled (regions optional)",
                 geojson_file)

    # ── 5. Write the unfiltered contract (X = raw counts) ────────────────
    Path(out_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)
    log.info("Wrote contract %s  (%d cells × %d genes, raw counts)",
             out_h5ad, adata.n_obs, adata.n_vars)

except Exception:
    log.error("FAILED for %s:\n%s", snakemake.params.sample_id, traceback.format_exc())
    raise
