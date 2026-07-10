#!/usr/bin/env python3
"""
Visium HD Public Dataset → AnnData (+ optional ROI H5ADs)
-----------------------------------------------------------
Input  : The 10 files 10x Genomics gives you when you download a Visium HD
         dataset (see REQUIRED_FILES below) — download manually, 10x blocks
         scripted downloads.
Output : One AnnData (.h5ad) with spatial coordinates for the chosen bin
         size, optionally split into per-ROI H5ADs using a GeoJSON.

Dependencies: numpy, scipy, h5py, pandas, anndata, geopandas, shapely
              pandas.read_parquet also needs pyarrow (or fastparquet):
                  pip install pyarrow
"""

import json
import shutil
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
import anndata as ad

# ── Paths ────────────────────────────────────────────────────────────────────
# Point this at the folder where you extracted/placed the 10 files 10x gives
# you (Visium_HD_Mouse_Brain_binned_outputs.tar.gz, _spatial.tar.gz, etc.)
DATA_DIR   = Path(r"D:\compressed_data\brain_visium")
SAMPLE     = "Visium_HD_Mouse_Brain"     # filename prefix used by 10x
OUTPUT_DIR = DATA_DIR.parent / "ROI_H5ADs"

BIN_SIZE      = "008um"       # one of "002um", "008um", "016um"
USE_FILTERED  = True          # filtered (tissue-covered) vs raw matrix

VISIUM_HD_DOWNLOAD_PAGE = (
    "https://www.10xgenomics.com/datasets"  # dataset-specific page varies
)

# The 10 files 10x Genomics ships for a Visium HD sample. Some are optional
# depending on whether you ran segmentation / generated a cell-level cloupe.
REQUIRED_FILES = [
    f"{SAMPLE}_binned_outputs.tar.gz",
    f"{SAMPLE}_spatial.tar.gz",
    f"{SAMPLE}_metrics_summary.csv",
]
OPTIONAL_FILES = [
    f"{SAMPLE}_barcode_mappings.parquet",
    f"{SAMPLE}_cloupe_008um.cloupe",
    f"{SAMPLE}_cloupe_cell.cloupe",
    f"{SAMPLE}_feature_slice.h5",
    f"{SAMPLE}_molecule_info.h5",
    f"{SAMPLE}_segmented_outputs.tar.gz",
    f"{SAMPLE}_web_summary.html",
]


# ════════════════════════════════════════════════════════════════════════════
# 0 ─ Check the data is present (10x Genomics blocks scripted downloads,
#     so the dataset must be fetched manually from the browser)
# ════════════════════════════════════════════════════════════════════════════

def ensure_data(path: Path) -> None:
    """
    Verify the Visium HD download folder exists and contains the expected
    files. 10x Genomics' download links sit behind bot protection that
    blocks urllib/requests, so this dataset can't be fetched
    programmatically — it must be downloaded manually from the dataset page.
    """
    if not path.exists() or not any(path.iterdir()):
        raise FileNotFoundError(
            f"[error] Data folder not found or empty: {path}\n\n"
            f"This dataset must be downloaded manually (10x Genomics blocks "
            f"scripted downloads):\n"
            f"  1. Open the dataset's page on {VISIUM_HD_DOWNLOAD_PAGE}\n"
            f"  2. Go to the 'Output and supplemental files' tab.\n"
            f"  3. Download the files listed there (binned_outputs.tar.gz, "
            f"spatial.tar.gz, etc.).\n"
            f"  4. Place them directly in:\n"
            f"       {path}\n"
        )

    missing_required = [f for f in REQUIRED_FILES if not (path / f).exists()]
    if missing_required:
        raise FileNotFoundError(
            f"[error] Data folder found but missing required file(s):\n"
            f"  " + "\n  ".join(missing_required) + "\n"
            f"  Folder checked: {path}\n"
            f"  Re-check your download from {VISIUM_HD_DOWNLOAD_PAGE} — "
            f"did all files finish downloading?\n"
        )

    missing_optional = [f for f in OPTIONAL_FILES if not (path / f).exists()]
    print(f"[info] Data folder found: {path}")
    if missing_optional:
        print(f"[info] Optional files not present (fine to skip): "
              f"{', '.join(missing_optional)}")


# ════════════════════════════════════════════════════════════════════════════
# 1 ─ Extract archives (only if not already extracted)
# ════════════════════════════════════════════════════════════════════════════

def _extract_tar(archive: Path, dest: Path) -> None:
    print(f"[step 1] Extracting {archive.name} → {dest} …")
    with tarfile.open(archive, "r:gz") as tf:
        tf.extractall(dest)


def ensure_extracted(data_dir: Path) -> tuple[Path, Path]:
    """
    Extracts *_binned_outputs.tar.gz and *_spatial.tar.gz next to
    themselves if not already extracted. Returns
    (binned_outputs_dir, spatial_dir).
    """
    binned_archive = data_dir / f"{SAMPLE}_binned_outputs.tar.gz"
    spatial_archive = data_dir / f"{SAMPLE}_spatial.tar.gz"

    binned_dir = data_dir / "binned_outputs"
    if not binned_dir.exists():
        _extract_tar(binned_archive, data_dir)
    else:
        print(f"[step 1] {binned_dir.name} already extracted, skipping.")

    spatial_dir = data_dir / "spatial"
    if not spatial_dir.exists():
        _extract_tar(spatial_archive, data_dir)
    else:
        print(f"[step 1] {spatial_dir.name} already extracted, skipping.")

    return binned_dir, spatial_dir


# ════════════════════════════════════════════════════════════════════════════
# 2 ─ Load a bin-resolution AnnData
# ════════════════════════════════════════════════════════════════════════════

def _read_10x_h5_matrix(h5_path: Path) -> ad.AnnData:
    """Read a Space Ranger filtered/raw_feature_bc_matrix.h5 into AnnData."""
    with h5py.File(h5_path, "r") as f:
        grp = f["matrix"]
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        shape = grp["shape"][:]          # (n_genes, n_barcodes)

        X = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()

        barcodes = grp["barcodes"][:].astype(str)
        feat = grp["features"]
        var = pd.DataFrame({
            "gene_ids": feat["id"][:].astype(str),
            "feature_types": feat["feature_type"][:].astype(str),
        }, index=feat["name"][:].astype(str))

    obs = pd.DataFrame(index=barcodes)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def load_visium_hd_anndata(binned_outputs_dir: Path,
                            bin_size: str = BIN_SIZE,
                            use_filtered: bool = USE_FILTERED) -> ad.AnnData:
    """
    Load one bin resolution (e.g. "008um") of a Visium HD run into AnnData,
    with spatial coordinates from tissue_positions.parquet attached to
    .obsm["spatial"] (pixel coordinates in the full-res image) and
    .obs["array_row"] / .obs["array_col"] for the bin grid position.
    """
    bin_dir = binned_outputs_dir / f"square_{bin_size}"
    if not bin_dir.exists():
        available = [p.name for p in binned_outputs_dir.glob("square_*")]
        raise FileNotFoundError(
            f"[error] No '{bin_dir.name}' folder under {binned_outputs_dir}. "
            f"Available resolutions: {available}"
        )

    matrix_name = "filtered_feature_bc_matrix.h5" if use_filtered \
        else "raw_feature_bc_matrix.h5"
    h5_path = bin_dir / matrix_name
    print(f"[step 2] Loading {matrix_name} ({bin_size}) …")
    adata = _read_10x_h5_matrix(h5_path)
    print(f"[step 2]   {adata.n_obs} bins × {adata.n_vars} genes")

    pos_path = bin_dir / "spatial" / "tissue_positions.parquet"
    try:
        positions = pd.read_parquet(pos_path)
    except ImportError as exc:
        raise ImportError(
            "reading tissue_positions.parquet requires pyarrow: "
            "pip install pyarrow"
        ) from exc

    positions = positions.set_index("barcode").loc[adata.obs_names]
    adata.obs["array_row"] = positions["array_row"].values
    adata.obs["array_col"] = positions["array_col"].values
    adata.obs["in_tissue"] = positions["in_tissue"].values.astype(bool)
    adata.obsm["spatial"] = positions[
        ["pxl_col_in_fullres", "pxl_row_in_fullres"]
    ].values.astype(float)

    scalefactors_path = bin_dir / "spatial" / "scalefactors_json.json"
    if scalefactors_path.exists():
        with open(scalefactors_path, "r") as fh:
            adata.uns["spatial_scalefactors"] = json.load(fh)

    adata.uns["bin_size_um"] = bin_size
    print(f"[step 2] AnnData: {adata.n_obs} bins × {adata.n_vars} genes "
          f"(bin_size={bin_size}, filtered={use_filtered})")
    return adata


# ════════════════════════════════════════════════════════════════════════════
# 3 ─ Optional: split by ROI (GeoJSON polygons drawn on the full-res image)
# ════════════════════════════════════════════════════════════════════════════

def find_geojson(data_dir: Path) -> Path | None:
    hits = list(data_dir.glob("*.geojson"))
    if not hits:
        print("[step 3] No GeoJSON found — skipping ROI split.")
        return None
    print(f"[step 3] Using GeoJSON: {hits[0].name}")
    return hits[0]


def split_by_roi(adata: ad.AnnData, geojson_path: Path) -> dict[str, ad.AnnData]:
    """
    Same approach as the Xenium script: point-in-polygon test between each
    bin's full-res pixel coordinate and each ROI polygon. Visium HD ROI
    annotations are normally drawn directly on the same full-res image the
    pixel coordinates come from, so — unlike Xenium — no pixel/micron
    rescaling is expected here. If you see 0 cells per ROI, check the
    printed bounds below for a scale/offset mismatch, same diagnostic idea
    as before.
    """
    import geopandas as gpd
    from shapely.geometry import Point

    gdf = gpd.read_file(str(geojson_path))

    from shapely.affinity import scale

    xy = adata.obsm["spatial"]

    gx_min, gy_min, gx_max, gy_max = gdf.total_bounds

    scale_factor = 1.0 / adata.uns["spatial_scalefactors"]["tissue_hires_scalef"]

    print(f"[step 3] Scaling GeoJSON by {scale_factor:.3f} (hires → fullres)")

    from shapely.affinity import scale

    gdf = gpd.read_file(str(geojson_path))

    print(f"[step 3] GeoJSON loaded: {len(gdf)} features")

    scale_factor = 1.0 / adata.uns["spatial_scalefactors"]["tissue_hires_scalef"]

    print(f"[step 3] Scaling GeoJSON by {scale_factor:.3f} (hires → fullres)")

    gdf["geometry"] = gdf.geometry.apply(
        lambda geom: scale(
            geom,
            xfact=scale_factor,
            yfact=scale_factor,
            origin=(0, 0),
        )
    )
    print(f"[step 3] GeoJSON loaded: {len(gdf)} features")

    xy = adata.obsm["spatial"]
    gx_min, gy_min, gx_max, gy_max = gdf.total_bounds
    print(f"[debug] GeoJSON bounds : x[{gx_min:.1f}, {gx_max:.1f}]  "
          f"y[{gy_min:.1f}, {gy_max:.1f}]")
    print(f"[debug] Bins bounds    : x[{xy[:,0].min():.1f}, {xy[:,0].max():.1f}]  "
          f"y[{xy[:,1].min():.1f}, {xy[:,1].max():.1f}]")

    points = gpd.GeoSeries([Point(x, y) for x, y in xy], crs=gdf.crs)

    roi_adatas: dict[str, ad.AnnData] = {}
    for i, row in gdf.iterrows():
        name = (row.get("name")
                 or (row.get("classification") or {}).get("name")
                 if isinstance(row.get("classification"), dict) else None)
        name = name or f"ROI_{i+1}"
        mask = points.within(row.geometry).values
        n = int(mask.sum())
        print(f"[step 3]   {name}: {n} bins")
        if n == 0:
            print(f"[warn]     No bins found in {name} — skipping.")
            continue
        roi_adatas[name] = adata[mask].copy()

    if not roi_adatas:
        raise RuntimeError(
            "No bins were assigned to any ROI. Check that the GeoJSON "
            "coordinates match the Visium HD full-res pixel space."
        )
    return roi_adatas


# ════════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════════

def main() -> None:
    ensure_data(DATA_DIR)
    binned_dir, spatial_dir = ensure_extracted(DATA_DIR)

    adata = load_visium_hd_anndata(binned_dir, BIN_SIZE, USE_FILTERED)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / f"{SAMPLE}_{BIN_SIZE}.h5ad"
    adata.write_h5ad(out_path)
    print(f"[done] Wrote {out_path}")

    geojson_path = find_geojson(DATA_DIR)
    if geojson_path is not None:
        roi_adatas = split_by_roi(adata, geojson_path)
        for name, roi_adata in roi_adatas.items():
            roi_path = OUTPUT_DIR / f"{SAMPLE}_{BIN_SIZE}_{name}.h5ad"
            roi_adata.write_h5ad(roi_path)
            print(f"[done] Wrote {roi_path}")


if __name__ == "__main__":
    main()
