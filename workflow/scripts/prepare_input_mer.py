"""
prepare_input_mer.py — MERSCOPE head: build the standardized UNFILTERED CONTRACT h5ad
from the Vizgen cell table + QuPath GeoJSON.
=====================================================================================
Analogous to prepare_input_x5k, but reads the two Vizgen CSVs directly (no zarr, no
squidpy) and reuses the coordinate metadata that generate_qupath_mer already computed.
Produces ONLY the contract:
  X                      raw integer counts (NO filtering / normalisation / clustering)
  obsm['spatial']        cell centroids, in microns FROM THE MOSAIC ORIGIN (image-origin
                         frame: raw center_x/y minus p0). A constant shift — invariant for
                         every downstream consumer (neighbourhoods, BANKSY, plotting) — that
                         aligns the cells with the embedded background and reduces the
                         GeoJSON round-trip to the Xenium head's bare-scalar form.
  uns['spatial'][sample] {images:{hires}, scalefactors}  (grey mosaic composite)
  obs                    sample, cell_id, region_annotation

QC, normalisation, PCA/UMAP and clustering are the CoreBlock's job — as for every head.
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
log = logging.getLogger("prepare_input_mer")

try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _here = os.getcwd()
sys.path.insert(0, _here)

try:
    import anndata as ad
    import pandas as pd
    from scipy.sparse import csr_matrix

    from merscope_io import find_merscope_files, read_cell_table

    merscope_dir = str(snakemake.params.merscope_dir)
    out_h5ad     = str(snakemake.output.h5ad)
    sample_id    = snakemake.params.sample_id
    geojson_dir  = snakemake.params.geojson_dir
    qupath_meta  = str(snakemake.input.qupath_meta)
    background   = str(snakemake.input.background)

    log.info("=" * 70)
    log.info("MERSCOPE contract for sample: %s", sample_id)
    log.info("=" * 70)

    # ── 1. Read the Vizgen cell table (RAW counts) ───────────────────────
    files = find_merscope_files(merscope_dir)
    counts, meta = read_cell_table(files["counts"], files["meta"])
    adata = ad.AnnData(
        X=csr_matrix(counts.values),
        obs=meta.copy(),
        var=pd.DataFrame(index=counts.columns.astype(str)),
    )
    adata.var_names_make_unique()
    adata.obs["sample"] = sample_id
    adata.obs["cell_id"] = adata.obs_names.astype(str)   # contract: cell_id == obs_names
    log.info("AnnData: %d cells × %d genes (raw, unfiltered)", adata.n_obs, adata.n_vars)

    # centroids in microns (Vizgen cell_metadata center_x / center_y)
    xcol = "center_x" if "center_x" in meta.columns else meta.columns[meta.columns.str.lower() == "center_x"][0]
    ycol = "center_y" if "center_y" in meta.columns else meta.columns[meta.columns.str.lower() == "center_y"][0]
    spatial = meta[[xcol, ycol]].to_numpy(dtype=float)

    # ── 2. Coordinate metadata (computed once by generate_qupath_mer) ─────
    with open(qupath_meta) as f:
        qm = json.load(f)
    p0 = np.asarray(qm["p0_micron"], dtype=float)
    scalef = float(qm["tissue_hires_scalef"])
    pixel_size_um = float(qm["pixel_size_um"])
    affine6 = list(qm["qupath_to_micron_affine"])

    # shift to the image-origin frame so pixel(0,0) <-> micron(0,0)
    adata.obsm["spatial"] = spatial - p0
    log.info("Shifted centroids to image-origin frame (p0 = %.2f, %.2f µm)", p0[0], p0[1])

    # ── 3. Embed the grey mosaic background (shares the QuPath grid) ──────
    grey = np.load(background)
    adata.uns["spatial"] = {
        sample_id: {
            "images": {"hires": grey},
            "scalefactors": {
                "tissue_hires_scalef": scalef,
                "spot_diameter_fullres": 10.0,
                "pixel_size_um": pixel_size_um,
            },
        }
    }
    adata.uns["morphology_channels"] = list(qm.get("channels", []))
    log.info("Embedded grey background %s (scalef=%.6f)", grey.shape, scalef)

    # ── 4. QuPath GeoJSON → region_annotation (optional) ─────────────────
    adata.obs["region_annotation"] = "Unlabeled"
    geojson_file = os.path.join(geojson_dir or "", f"{sample_id}_morphology.geojson")
    if geojson_dir and os.path.isfile(geojson_file):
        log.info("Reading QuPath GeoJSON: %s", geojson_file)
        try:
            import geopandas as gpd
            annot_gdf = gpd.read_file(geojson_file)
            annot_gdf["region"] = annot_gdf["classification"].apply(
                lambda x: json.loads(x).get("name", "Unlabeled")
                if isinstance(x, str) else "Unlabeled")
            # QuPath pixel polygons -> image-origin micron (affine has no offset; for the
            # usual diagonal transform this equals x5k's geometry.scale(px_to_um, origin=0))
            annot_gdf["geometry"] = annot_gdf["geometry"].affine_transform(affine6)
            sp = adata.obsm["spatial"]
            cell_gdf = gpd.GeoDataFrame(
                {"cell_id": adata.obs["cell_id"].values},
                geometry=gpd.points_from_xy(sp[:, 0], sp[:, 1]))
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
    log.error("FAILED for %s:\n%s",
              getattr(snakemake.params, "sample_id", "?"), traceback.format_exc())
    raise
