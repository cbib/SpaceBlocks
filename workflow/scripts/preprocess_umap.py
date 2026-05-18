"""
preprocess_umap.py – Load, QC, normalise, embed
=================================================
Discovers Space Ranger output files from the SR output directory.
Uses tissue_positions.parquet + barcode_mappings.parquet for spatial
coordinates (post-rotation, aligned with tissue_hires_image.png).
Uses the **raw** cell matrix.
"""

import gc
import json
import logging
import os
import shutil
import sys
import traceback
from pathlib import Path
from glob import glob
from datetime import datetime

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import scanpy as sc
import scanpy.external as sce


# ── Logging ──────────────────────────────────────────────────────────────────
log_handlers = [logging.StreamHandler(sys.stderr)]
if hasattr(snakemake, "log"):
    if snakemake.log.out:
        Path(snakemake.log.out).parent.mkdir(parents=True, exist_ok=True)
        log_handlers.append(logging.FileHandler(snakemake.log.out, mode="w"))
    if snakemake.log.err:
        Path(snakemake.log.err).parent.mkdir(parents=True, exist_ok=True)
        sys.stderr = open(snakemake.log.err, "w")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=log_handlers,
)
log = logging.getLogger("preprocess_umap")

# ── Parameters ───────────────────────────────────────────────────────────────
sample_id      = snakemake.params.sample_id
sr_outdir      = snakemake.params.sr_outdir
geojson_path   = snakemake.params.geojson_path
MIN_COUNTS     = int(snakemake.params.min_counts)
MIN_CELLS      = int(snakemake.params.min_cells)
MIN_GENES      = int(snakemake.params.min_genes)
N_NEIGHBORS    = int(snakemake.params.n_neighbors)
N_PCS          = int(snakemake.params.n_pcs)
RES_SCAN_MIN   = float(snakemake.params.resolution_scan_min)
RES_SCAN_MAX   = float(snakemake.params.resolution_scan_max)
RES_SCAN_STEP  = float(snakemake.params.resolution_scan_step)
RANDOM_SEED    = int(snakemake.params.random_seed)
USE_PRECOMPUTED = bool(snakemake.params.use_precomputed)

out_adata      = str(snakemake.output.adata)
out_metadata   = str(snakemake.output.metadata)
out_qupath_img = str(snakemake.output.qupath_image)
output_dir     = str(Path(out_adata).parent)


# ── File discovery helper ────────────────────────────────────────────────────
def find_sr_file(sr_outdir, *patterns):
    """Search for a file under the Space Ranger output directory.
    Tries each glob pattern in order and returns the first match."""
    for pattern in patterns:
        matches = glob(os.path.join(sr_outdir, pattern))
        if matches:
            return matches[0]
    tried = ", ".join(patterns)
    raise FileNotFoundError(
        f"Could not find any of [{tried}] under {sr_outdir}"
    )


def _resolution_range(rmin, rmax, step):
    n = round((rmax - rmin) / step)
    return [round(rmin + i * step, 10) for i in range(n + 1)]


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
    log.info("Preprocessing sample: %s", sample_id)
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

    hires_image_path = find_sr_file(
        sr_outdir,
        "outs/segmented_outputs/spatial/tissue_hires_image.png",
        "outs/spatial/tissue_hires_image.png",
        "segmented_outputs/spatial/tissue_hires_image.png",
        "spatial/tissue_hires_image.png",
    )

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

    del tissue_pos, barcode_map, cell_centroids
    gc.collect()

    # ── 4. QC ────────────────────────────────────────────────────────────
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hb"], inplace=True, log1p=False)

    log.info("QC: min_counts=%d, min_cells=%d, min_genes=%d",
             MIN_COUNTS, MIN_CELLS, MIN_GENES)
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)

    log.info("After QC: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # ── 5. QuPath region annotations (optional) ─────────────────────────
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
        region_labels = region_labels.replace("Necrosis", "Bubble")
        adata.obs["region_annotation"] = region_labels

        sc.pl.spatial(adata, color="region_annotation", spot_size=20,
                      title="Annotated Regions", library_id=sample_id)
        plt.savefig(os.path.join(output_dir, f"region_annotation_{sample_id}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
    else:
        log.warning("No QuPath annotation found for '%s'.", sample_id)
        adata.obs["region_annotation"] = "Unlabeled"

    # ── 6. Normalise ─────────────────────────────────────────────────────
    log.info("Normalising (target_sum=None + log1p) …")
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)

    # ── 7. PCA → UMAP ─────────────────────────────────────────
    log.info("PCA (all genes) → UMAP (seed=%d) …", RANDOM_SEED)
    sc.pp.pca(adata, use_highly_variable=False, random_state=RANDOM_SEED)
    adata.obsm["X_pca_original"] = adata.obsm["X_pca"].copy()
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)

    # ── 8. Leiden at all resolutions ─────────────────────────────────────
    resolutions = _resolution_range(RES_SCAN_MIN, RES_SCAN_MAX, RES_SCAN_STEP)
    leiden_keys = [f"leiden_{str(r).replace('.', '_')}" for r in resolutions]

    if USE_PRECOMPUTED and os.path.isfile(out_metadata):
        log.info("Reloading precomputed clusters from %s …", out_metadata)
        saved_meta = pd.read_csv(out_metadata, sep="\t", index_col=0, comment="#")
        for key in leiden_keys:
            if key in saved_meta.columns:
                adata.obs[key] = saved_meta[key].reindex(adata.obs_names).astype(str).astype("category")
                log.info("  Loaded %s: %d clusters", key, adata.obs[key].nunique())
            else:
                log.warning("  %s not in saved metadata, computing …", key)
                res = float(key.replace("leiden_", "").replace("_", "."))
                sc.tl.leiden(adata, resolution=res, key_added=key, random_state=RANDOM_SEED)
    else:
        log.info("Computing Leiden for %d resolutions: %s", len(resolutions), list(resolutions))
        for res in resolutions:
            key = f"leiden_{str(res).replace('.', '_')}"
            sc.tl.leiden(adata, resolution=res, key_added=key, random_state=RANDOM_SEED)
            log.info("  res=%.1f → %d clusters", res, adata.obs[key].nunique())

    # ── 9. Save ──────────────────────────────────────────────────────────
    log.info("Saving adata → %s", out_adata)
    Path(out_adata).parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_adata)

    # Save metadata TSV with cluster assignments for reproducibility
    log.info("Saving metadata TSV → %s", out_metadata)
    meta_cols = ["sample", "region_annotation"] + [k for k in leiden_keys if k in adata.obs.columns]
    meta_df = adata.obs[meta_cols].copy()
    Path(out_metadata).parent.mkdir(parents=True, exist_ok=True)
    with open(out_metadata, "w") as f:
        f.write(f"# Generated: {datetime.now().isoformat()}\n")
        f.write(f"# Sample: {sample_id}\n")
        f.write(f"# Scanpy: {sc.__version__}\n")
        f.write(f"# Random seed: {RANDOM_SEED}\n")
        f.write(f"# n_cells: {adata.n_obs}\n")
        meta_df.to_csv(f, sep="\t")

    log.info("Copying hires image → %s", out_qupath_img)
    Path(out_qupath_img).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(hires_image_path, out_qupath_img)

    log.info("Preprocessing complete for %s.", sample_id)

except Exception:
    log.error("FAILED for %s:\n%s", sample_id, traceback.format_exc())
    raise
