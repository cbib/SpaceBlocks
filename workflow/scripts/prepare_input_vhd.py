"""
prepare_input_vhd.py — Visium HD HEAD (step 3).
===========================================
Builds the UNFILTERED contract h5ad that the common core consumes. All
technology-specific work lives here:
  * discover Space Ranger output files (tree varies across versions),
  * load the **raw** cell matrix,
  * compute cell centroids from tissue_positions.parquet + barcode_mappings.parquet
    (post-rotation, aligned with tissue_hires_image.png),
  * embed the hires image + scalefactors in uns["spatial"] for sc.pl.spatial,
  * join the QuPath geojson into obs["region_annotation"].

NO filtering and NO normalization — X stays raw counts. Output contract:
X = raw counts; obsm["spatial"] = (x, y); uns["spatial"] = image + scalefactors;
obs has "sample", "cell_id", "region_annotation".
"""
import gc
import json
import logging
import os
import sys
import traceback
from glob import glob
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from PIL import Image
import scanpy as sc

# ── Logging ──────────────────────────────────────────────────────────────────
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
                    datefmt="%Y-%m-%d %H:%M:%S", handlers=log_handlers)
log = logging.getLogger("prepare_input")

# ── Parameters ───────────────────────────────────────────────────────────────
sample_id    = snakemake.params.sample_id
sr_outdir    = snakemake.params.sr_outdir
geojson_path = snakemake.params.geojson_path
out_h5ad     = str(snakemake.output.h5ad)


# ── File discovery helper ────────────────────────────────────────────────────
# Space Ranger file discovery, shared with the other VHD head step (no drift).
try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:                      # very old Snakemake
    _here = os.getcwd()
sys.path.insert(0, _here)
from vhd_sr import find_sr_file


def _normalise_cell_id(raw_id):
    """
    Normalise cell IDs to a canonical integer string for matching.

    Handles all known formats from Space Ranger:
      'cellid_000000123-1' → '123'    (adata barcodes)
      'cellid_000000123'   → '123'    (barcode_mappings)
      '123'                → '123'    (plain numeric)
      'cellid_000000000'   → '0'      (unassigned — will be filtered)
    """
    s = str(raw_id)
    # Strip trailing -1 suffix (only the last one)
    if s.endswith("-1"):
        s = s[:-2]
    # Strip cellid_ prefix
    if s.startswith("cellid_"):
        s = s[7:]
    # Convert to int then back to string to remove leading zeros
    try:
        return str(int(s))
    except ValueError:
        return s


try:
    log.info("=" * 70)
    log.info("Preparing input contract: %s", sample_id)
    log.info("SR output dir: %s", sr_outdir)
    log.info("=" * 70)

    # ── 1. Discover Space Ranger output files ────────────────────────────
    cell_matrix_path = find_sr_file(
        sr_outdir,
        "outs/segmented_outputs/raw_feature_cell_matrix.h5",
        "segmented_outputs/raw_feature_cell_matrix.h5",
        "outs/raw_feature_cell_matrix.h5",
    )

    tissue_pos_path = find_sr_file(
        sr_outdir,
        "outs/binned_outputs/square_002um/spatial/tissue_positions.parquet",
        "binned_outputs/square_002um/spatial/tissue_positions.parquet",
    )

    barcode_map_path = find_sr_file(
        sr_outdir,
        "outs/barcode_mappings.parquet",
        "barcode_mappings.parquet",
    )

    # canonical copy declared as an input (generate_qupath_vhd already discovered it) —
    # no second glob of the SR tree, so the two head steps can't diverge on the image.
    hires_image_path = str(snakemake.input.hires_png)

    scalefactors_path = find_sr_file(
        sr_outdir,
        "outs/segmented_outputs/spatial/scalefactors_json.json",
        "outs/spatial/scalefactors_json.json",
        "segmented_outputs/spatial/scalefactors_json.json",
        "spatial/scalefactors_json.json",
    )

    log.info("Discovered files:")
    log.info("  cell_matrix:      %s", cell_matrix_path)
    log.info("  tissue_positions: %s", tissue_pos_path)
    log.info("  barcode_mappings: %s", barcode_map_path)
    log.info("  hires_image:      %s", hires_image_path)
    log.info("  scalefactors:     %s", scalefactors_path)

    # ── 2. Load cell matrix ──────────────────────────────────────────────
    log.info("Loading raw cell matrix …")
    adata = sc.read_10x_h5(cell_matrix_path)
    adata.var_names_make_unique()
    adata.obs["sample"] = sample_id
    log.info("  Loaded %d cells, %d genes", adata.n_obs, adata.n_vars)

    # Log barcode format examples
    log.info("  Example adata barcodes: %s", list(adata.obs_names[:3]))

    with open(scalefactors_path) as f:
        scale_data = json.load(f)
    hires_scale = scale_data["tissue_hires_scalef"]
    log.info("  hires_scalef = %.6f", hires_scale)

    hires_img = np.array(Image.open(hires_image_path))

    # ── 3. Spatial coordinates from tissue_positions + barcode_mappings ──
    log.info("Computing cell centroids from tissue_positions + barcode_mappings …")
    tissue_pos = pd.read_parquet(tissue_pos_path)
    tissue_pos = tissue_pos.set_index("barcode")

    barcode_map = pd.read_parquet(barcode_map_path)
    log.info("  barcode_mappings columns: %s", list(barcode_map.columns))

    # Find cell_id column (case-insensitive)
    cell_id_col = None
    for col in barcode_map.columns:
        if col.lower() == "cell_id":
            cell_id_col = col
            break
    if cell_id_col is None:
        raise KeyError(f"No cell_id column in barcode_mappings. Columns: {list(barcode_map.columns)}")

    # Find square_002um column (case-insensitive)
    bin_col = None
    for col in barcode_map.columns:
        if col.lower() == "square_002um":
            bin_col = col
            break
    if bin_col is None:
        raise KeyError(f"No square_002um column in barcode_mappings. Columns: {list(barcode_map.columns)}")

    log.info("  Using columns: cell_id='%s', bin='%s'", cell_id_col, bin_col)
    log.info("  Example cell_ids: %s", list(barcode_map[cell_id_col].dropna().head(5)))

    # Normalise cell_ids and filter out unassigned bins (cell_id = 0)
    barcode_map["_norm_cell_id"] = barcode_map[cell_id_col].astype(str).map(_normalise_cell_id)

    # Filter: remove bins not assigned to a real cell (id = "0" or NaN)
    barcode_map = barcode_map[
        barcode_map["_norm_cell_id"].notna()
        & (barcode_map["_norm_cell_id"] != "0")
        & (barcode_map["_norm_cell_id"] != "nan")
    ].copy()
    log.info("  After filtering unassigned bins: %d bins", len(barcode_map))

    # Join bin positions with cell assignments
    barcode_map = barcode_map.merge(
        tissue_pos[["pxl_row_in_fullres", "pxl_col_in_fullres"]],
        left_on=bin_col,
        right_index=True,
        how="inner",
    )

    # Centroid = mean position of constituent 2µm bins per cell
    cell_centroids = (
        barcode_map
        .groupby("_norm_cell_id")[["pxl_row_in_fullres", "pxl_col_in_fullres"]]
        .mean()
    )
    log.info("  Computed centroids for %d cells from %d bins",
             len(cell_centroids), len(barcode_map))
    log.info("  Example centroid IDs: %s", list(cell_centroids.index[:5]))

    # Normalise adata barcodes the same way
    adata.obs["_norm_cell_id"] = adata.obs_names.map(_normalise_cell_id)
    log.info("  Example normalised adata IDs: %s", list(adata.obs["_norm_cell_id"][:5]))

    # Match
    has_coords = adata.obs["_norm_cell_id"].isin(cell_centroids.index)
    n_matched = has_coords.sum()
    n_dropped = (~has_coords).sum()
    log.info("  Matched: %d, Unmatched: %d", n_matched, n_dropped)

    if n_matched == 0:
        # Debug: show what the IDs look like
        log.error("  ZERO matches! Dumping ID samples for debugging:")
        log.error("  adata normalised IDs (first 10): %s",
                  list(adata.obs["_norm_cell_id"].head(10)))
        log.error("  centroid index (first 10): %s",
                  list(cell_centroids.index[:10]))
        raise RuntimeError(
            f"No cells matched between adata barcodes and barcode_mappings. "
            f"adata IDs look like: {list(adata.obs['_norm_cell_id'].head(3))}, "
            f"centroid IDs look like: {list(cell_centroids.index[:3])}"
        )

    if n_dropped > 0:
        log.warning("  Dropping %d cells without spatial coordinates (%.1f%%)",
                     n_dropped, 100 * n_dropped / adata.n_obs)
    adata = adata[has_coords].copy()

    # Store spatial coords: scanpy expects (x, y) = (pxl_col, pxl_row)
    matched = cell_centroids.loc[adata.obs["_norm_cell_id"].values]
    adata.obsm["spatial"] = np.column_stack([
        matched["pxl_col_in_fullres"].values,
        matched["pxl_row_in_fullres"].values,
    ])

    # Store image + scalefactors in uns for sc.pl.spatial
    adata.uns["spatial"] = {
        sample_id: {
            "images": {"hires": hires_img},
            "scalefactors": {
                "tissue_hires_scalef": hires_scale,
                "spot_diameter_fullres": scale_data.get("spot_diameter_fullres", 1.0),
            },
        }
    }

    # Prefix obs index with sample name
    adata.obs.index = sample_id + "_" + adata.obs.index.astype(str)
    adata.obs["cell_id"] = adata.obs.index
    log.info("  %d cells with registered spatial coordinates", adata.n_obs)

    # Drop the helper column; keep the object lean before the geojson join
    del adata.obs["_norm_cell_id"]
    del tissue_pos, barcode_map, cell_centroids
    gc.collect()

    # ── 4. QuPath region annotations (optional) ─────────────────────────
    qupath_geojson = os.path.join(
        geojson_path, f"{sample_id}_tissue_hires_image.geojson"
    )
    has_annotations = os.path.isfile(qupath_geojson)

    if has_annotations:
        log.info("QuPath annotation found: %s", qupath_geojson)
        annot_gdf = gpd.read_file(qupath_geojson)

        def _get_label(x):
            try:
                return json.loads(x)["name"]
            except (TypeError, json.JSONDecodeError):
                return "Unlabeled"

        annot_gdf["region"] = annot_gdf["classification"].apply(_get_label)

        spatial = adata.obsm["spatial"]
        coords_df = pd.DataFrame(
            {"x": spatial[:, 0], "y": spatial[:, 1]}, index=adata.obs_names
        )
        coords_df["cell_id"] = adata.obs_names.values
        # geojson is drawn on the hires image → scale full-res coords to hires
        coords_df["x_hires"] = coords_df["x"] * hires_scale
        coords_df["y_hires"] = coords_df["y"] * hires_scale

        cell_gdf = gpd.GeoDataFrame(
            coords_df,
            geometry=gpd.points_from_xy(coords_df["x_hires"], coords_df["y_hires"]),
            crs=annot_gdf.crs,
        )
        joined = gpd.sjoin(cell_gdf, annot_gdf[["region", "geometry"]],
                           how="left", predicate="within")
        joined = joined.drop_duplicates(subset="cell_id", keep="first")
        region_labels = joined.set_index("cell_id")["region"]
        region_labels = region_labels.reindex(adata.obs_names, fill_value="Unlabeled")
        # Necrosis polygons are air bubbles / artefacts → rename for downstream
        region_labels = region_labels.replace("Necrosis", "Bubble")
        adata.obs["region_annotation"] = region_labels
        log.info("  region_annotation value counts:\n%s",
                 adata.obs["region_annotation"].value_counts().to_string())
    else:
        log.warning("No QuPath annotation found for '%s'. "
                    "Setting region_annotation = 'Unlabeled'.", sample_id)
        adata.obs["region_annotation"] = "Unlabeled"

    # ── 5. Save UNFILTERED contract h5ad ─────────────────────────────────
    Path(out_h5ad).parent.mkdir(parents=True, exist_ok=True)
    log.info("Saving unfiltered contract h5ad → %s", out_h5ad)
    log.info("  %d cells, %d genes, X = raw counts (no filtering, no normalization)",
             adata.n_obs, adata.n_vars)
    adata.write(out_h5ad)

    log.info("prepare_input complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
