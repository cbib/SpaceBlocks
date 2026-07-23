#!/usr/bin/env python3
"""
Xenium 5k Public Dataset → ROI H5ADs + TSV metadata
------------------------------------------------------
Input  : Xenium Prime Human Skin FFPE output folder (zarr)
Output : One H5AD per GeoJSON ROI  +  metadata TSV
"""

import os
import sys
import json
from pathlib import Path

import numpy as np
import geopandas as gpd
import anndata as ad
from shapely.geometry import Point
from shapely.affinity import scale as shapely_affinity_scale

# ── 10x Genomics public dataset (download manually — see ensure_data below) ─
XENIUM_DOWNLOAD_PAGE = (
    "https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-skin"
)

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR   = Path(r"D:\Xenium designer\Xenium_Prime_Human_Skin_FFPE_xe_outs")
OUTPUT_DIR = DATA_DIR.parent / "ROI_H5ADs"          # one level up
SLIDE_ID   = "slide_1"
AREA_ID    = "A1"

# Files that must be present for the loader to work. Xenium bundles ship
# cells.zarr either as a folder or as cells.zarr.zip (both are handled by
# the loader below) — accept either form here. gene_panel.json is
# similarly sometimes named genes.json depending on export version.
REQUIRED_XENIUM_FILE_GROUPS = [
    ["cells.zarr", "cells.zarr.zip"],
    ["gene_panel.json", "genes.json"],
]


# ════════════════════════════════════════════════════════════════════════════
# 0 ─ Check the data is present (10x Genomics blocks scripted downloads,
#     so the dataset must be fetched manually from the browser)
# ════════════════════════════════════════════════════════════════════════════

def ensure_data(path: Path) -> None:
    """
    Verify the Xenium output folder exists and contains the expected files.
    10x Genomics' download links (cf.10xgenomics.com) sit behind bot
    protection that blocks urllib/requests, so this dataset can't be
    fetched programmatically — it must be downloaded manually.
    """
    if not path.exists() or not any(path.iterdir()):
        raise FileNotFoundError(
            f"[error] Data folder not found or empty: {path}\n\n"
            f"This dataset must be downloaded manually (10x Genomics blocks "
            f"scripted downloads):\n"
            f"  1. Open: {XENIUM_DOWNLOAD_PAGE}\n"
            f"  2. Go to the 'Output and supplemental files' tab.\n"
            f"  3. Download the Xenium Output Bundle (.zip) and extract it.\n"
            f"  4. Point DATA_DIR at the extracted folder (the one "
            f"containing cells.zarr, gene_panel.json, etc.), or move/rename "
            f"the extracted folder to:\n"
            f"       {path}\n"
        )

    missing_groups = [
        group[0] for group in REQUIRED_XENIUM_FILE_GROUPS
        if not any((path / name).exists() for name in group)
    ]
    if missing_groups:
        raise FileNotFoundError(
            f"[error] Data folder found but missing expected file(s): "
            f"{', '.join(missing_groups)}\n"
            f"  Folder checked: {path}\n"
            f"  Make sure you extracted the full Xenium Output Bundle from:\n"
            f"    {XENIUM_DOWNLOAD_PAGE}\n"
        )

    print(f"[info] Data folder found: {path}")


# ════════════════════════════════════════════════════════════════════════════
# 1 ─ Build AnnData from Xenium zarr output
# ════════════════════════════════════════════════════════════════════════════

def load_xenium_anndata(data_dir: Path) -> ad.AnnData:
    """
    Load Xenium output into AnnData.
    Uses spatialdata_io when available, otherwise falls back to
    manual zarr / parquet loading so the script is self-contained.
    """
    try:
        import spatialdata_io
        print("[step 1] Loading via spatialdata_io …")
        sdata = spatialdata_io.xenium(str(data_dir))

        # Pull the cell-by-gene table (the default table in Xenium sdatas)
        table_key = [k for k in sdata.tables][0]
        adata = sdata.tables[table_key]

        # Attach (x, y) spatial coordinates to obsm if not already there
        if "spatial" not in adata.obsm:
            cells = sdata["cell_circles"]           # geopandas GeoDataFrame
            xy = np.column_stack(
                [cells.geometry.x.values, cells.geometry.y.values]
            )
            adata.obsm["spatial"] = xy.astype(np.float32)

        print(f"[step 1] AnnData: {adata.n_obs} cells × {adata.n_vars} genes")
        return adata

    except Exception as e_spio:
        print(f"[step 1] spatialdata_io not available ({e_spio}); "
              "falling back to manual loader …")
        return _manual_xenium_loader(data_dir)


def _manual_xenium_loader(data_dir: Path) -> ad.AnnData:
    import zarr
    import pandas as pd
    import scipy.sparse as sp

    # ── 1. Cell metadata — ALL cells from cells.zarr ─────────────────────
    cells_zarr = data_dir / "cells.zarr.zip"
    store = (zarr.ZipStore(str(cells_zarr), mode="r")
             if cells_zarr.exists()
             else zarr.DirectoryStore(str(data_dir / "cells.zarr")))
    root = zarr.open(store, mode="r")

    cell_id_raw  = np.array(root["cell_id"])
    all_cell_ids = (cell_id_raw[:, 0] if cell_id_raw.ndim == 2
                    else cell_id_raw).astype(str)

    cell_summary = np.array(root["cell_summary"])
    cs_attrs     = dict(root["cell_summary"].attrs)
    col_names    = cs_attrs.get("column_names") or [
        "cell_centroid_x", "cell_centroid_y", "cell_area",
        "nucleus_centroid_x", "nucleus_centroid_y", "nucleus_area",
        "z_level", "nucleus_count"
    ][:cell_summary.shape[1]]

    def _find_col(names, *candidates):
        for c in candidates:
            if c in names:
                return names.index(c)
        raise KeyError(f"Cannot find any of {candidates} in {names}")

    x_col = _find_col(col_names, "cell_centroid_x", "x_centroid", "x")
    y_col = _find_col(col_names, "cell_centroid_y", "y_centroid", "y")

    # Full obs for ALL 112551 cells
    obs_full = pd.DataFrame(cell_summary, columns=col_names, index=all_cell_ids)
    obs_full.insert(0, "cell_id", all_cell_ids)
    obs_full["x_centroid"] = cell_summary[:, x_col].astype(np.float32)
    obs_full["y_centroid"] = cell_summary[:, y_col].astype(np.float32)
    print(f"[step 1]   {len(all_cell_ids)} total cells in cells.zarr")

    # ── 2. Feature names from gene_panel.json ────────────────────────────
    gene_panel_path = data_dir / "gene_panel.json"
    if not gene_panel_path.exists():
        raise FileNotFoundError(f"gene_panel.json not found in {data_dir}")
    with open(gene_panel_path, "r", encoding="utf-8") as fh:
        gene_panel = json.load(fh)
    targets  = gene_panel.get("payload", gene_panel).get("targets", [])
    features = [t.get("gene_name") or t.get("gene_identifier", f"gene_{i}")
                for i, t in enumerate(targets)]
    print(f"[step 1]   {len(features)} features from gene_panel.json")

    # ── 3. Expression matrix (10018 cells with QC-passing expression) ────
    mat_candidates = [
        data_dir / "cell_feature_matrix.zarr",
        data_dir / "cell_feature_matrix.zarr.zip",
    ]
    mat_path = next((p for p in mat_candidates if p.exists()), None)
    if mat_path is None:
        raise FileNotFoundError(
            f"Cannot find cell_feature_matrix.zarr[.zip] in {data_dir}"
        )

    mat_store = (zarr.ZipStore(str(mat_path), mode="r")
                 if str(mat_path).endswith(".zip")
                 else zarr.DirectoryStore(str(mat_path)))
    mat_root = zarr.open(mat_store, mode="r")
    cf       = mat_root["cell_features"]

    data_arr    = np.array(cf["data"],    dtype=np.float32)
    indices_arr = np.array(cf["indices"], dtype=np.int32)
    indptr_arr  = np.array(cf["indptr"],  dtype=np.int32)

    n_cells_mat  = len(indptr_arr) - 1
    mat_id_raw   = np.array(cf["cell_id"])
    mat_cell_ids = (mat_id_raw[:, 0] if mat_id_raw.ndim == 2
                    else mat_id_raw)[:n_cells_mat].astype(str)
    print(f"[step 1]   {n_cells_mat} cells with expression data (QC-passing)")

    X_mat = sp.csr_matrix(
        (data_arr, indices_arr, indptr_arr),
        shape=(n_cells_mat, len(features)),
    )

    # ── 4. Build unified AnnData over ALL cells ───────────────────────────
    # Strategy: 
    #   - reorder obs_full so expression cells come first (matching X_mat row order)
    #   - append zero rows for non-expression cells
    #   - then reindex back to original order

    n_all   = len(all_cell_ids)
    n_genes = len(features)

    # Which cells have expression, in matrix row order
    mat_id_to_row = pd.Series(np.arange(n_cells_mat), index=mat_cell_ids)

    # Boolean mask over all cells
    has_expr = obs_full.index.isin(set(mat_cell_ids))
    n_expr   = has_expr.sum()
    n_noexpr = n_all - n_expr
    print(f"[step 1]   {n_expr} cells with expression, "
          f"{n_noexpr} cells without (zeros)")

    # Reorder obs: expression cells first (in matrix row order), then the rest
    obs_expr   = obs_full.loc[mat_cell_ids]                    # n_cells_mat rows, matrix order
    obs_noexpr = obs_full.loc[~has_expr]                       # remaining cells

    # Corresponding X blocks
    # X_mat rows already align with obs_expr
    X_zeros = sp.csr_matrix((n_noexpr, n_genes), dtype=np.float32)
    X_stacked = sp.vstack([X_mat, X_zeros], format="csr")      # (n_all, n_genes)

    obs_stacked = pd.concat([obs_expr, obs_noexpr])

    # Add expression flag
    obs_stacked["has_expression"] = (
        [True]  * n_expr +
        [False] * n_noexpr
    )

    # Restore original cell order (optional but clean)
    original_order = list(all_cell_ids)
    obs_stacked    = obs_stacked.reindex(original_order)

    # Build a reindex map to reorder X_stacked rows to match obs_stacked
    stacked_index  = list(obs_expr.index) + list(obs_noexpr.index)
    pos_map        = {cid: i for i, cid in enumerate(stacked_index)}
    new_row_order  = np.array([pos_map[cid] for cid in original_order])
    X_full         = X_stacked[new_row_order, :]

    print(f"[step 1]   X_full shape: {X_full.shape}")

    var = pd.DataFrame(index=pd.Index(features, name="gene"))

    adata = ad.AnnData(X=X_full, obs=obs_stacked, var=var)
    adata.obsm["spatial"] = np.column_stack(
        [obs_stacked["x_centroid"].values,
         obs_stacked["y_centroid"].values]
    ).astype(np.float32)

    print(f"[step 1] AnnData: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata
# ════════════════════════════════════════════════════════════════════════════
# 2 ─ Ingest GeoJSON and split AnnData into ROIs
# ════════════════════════════════════════════════════════════════════════════

def find_geojson(data_dir: Path) -> Path:
    """Return the first .geojson file found in data_dir."""
    hits = list(data_dir.glob("*.geojson"))
    if not hits:
        raise FileNotFoundError(
            f"No .geojson file found in {data_dir}. "
            "Export your QuPath annotations as GeoJSON and place them there."
        )
    if len(hits) > 1:
        print(f"[warn] Multiple GeoJSON files found; using: {hits[0].name}")
    return hits[0]


def _detect_pixel_size(data_dir: Path, default: float = 0.2125) -> float:
    """
    Try to read the micron-per-pixel scale factor from Xenium's own
    metadata (experiment.xenium is a JSON file despite the extension).
    Falls back to the standard Xenium value if not found.
    """
    candidates = [
        data_dir / "experiment.xenium",
        data_dir / "experiment.xenium.json",
    ]
    for path in candidates:
        if path.exists():
            try:
                with open(path, "r", encoding="utf-8") as fh:
                    meta = json.load(fh)
                for key in ("pixel_size", "pixelSize", "micronsPerPixel"):
                    if key in meta:
                        val = float(meta[key])
                        print(f"[step 2]   Found pixel size in {path.name}: "
                              f"{val} µm/px")
                        return val
            except Exception:
                pass
    print(f"[step 2]   No pixel size metadata found; assuming default "
          f"{default} µm/px")
    return default


def split_by_roi(adata: ad.AnnData,
                 geojson_path: Path,
                 data_dir: Path | None = None) -> dict[str, ad.AnnData]:
    """
    For each polygon feature in the GeoJSON, collect all cells whose
    (x_centroid, y_centroid) falls inside, and return a sub-AnnData.

    The ROI name comes from (in priority order):
        feature["properties"]["name"]
        feature["properties"]["classification"]["name"]
        feature["properties"]["objectType"]
        "ROI_{i+1}"
    """
    gdf = gpd.read_file(str(geojson_path))
    print(f"[step 2] GeoJSON loaded: {len(gdf)} features")

    xy = adata.obsm["spatial"]          # shape (n_cells, 2)

    # ── Detect / correct coordinate-space mismatches (pixels vs microns) ───
    gx_min, gy_min, gx_max, gy_max = gdf.total_bounds
    cx_min, cy_min = xy[:, 0].min(), xy[:, 1].min()
    cx_max, cy_max = xy[:, 0].max(), xy[:, 1].max()
    print(f"[debug] GeoJSON bounds : x[{gx_min:.1f}, {gx_max:.1f}]  "
          f"y[{gy_min:.1f}, {gy_max:.1f}]")
    print(f"[debug] Cells bounds   : x[{cx_min:.1f}, {cx_max:.1f}]  "
          f"y[{cy_min:.1f}, {cy_max:.1f}]")

    g_span_x = gx_max - gx_min
    c_span_x = cx_max - cx_min
    ratio = g_span_x / c_span_x if c_span_x else 1.0

    # If the GeoJSON footprint is much larger than the cells' footprint,
    # assume it's in pixel space and rescale it into microns.
    if ratio > 1.5 or ratio < (1 / 1.5):
        pixel_size = _detect_pixel_size(data_dir) if data_dir else 0.2125
        print(f"[step 2]   Coordinate scale mismatch detected "
              f"(ratio={ratio:.2f}); rescaling GeoJSON by {pixel_size} "
              f"(pixels → microns)")
        gdf["geometry"] = gdf["geometry"].apply(
            lambda geom: shapely_affinity_scale(
                geom, xfact=pixel_size, yfact=pixel_size, origin=(0, 0)
            )
        )
        gx_min, gy_min, gx_max, gy_max = gdf.total_bounds
        print(f"[debug] GeoJSON bounds (rescaled): "
              f"x[{gx_min:.1f}, {gx_max:.1f}]  y[{gy_min:.1f}, {gy_max:.1f}]")
    # ──────────────────────────────────────────────────────────────────────
    points = gpd.GeoSeries(
        [Point(x, y) for x, y in xy],
        crs=gdf.crs                     # match CRS (usually None for Xenium)
    )

    roi_adatas: dict[str, ad.AnnData] = {}

    for i, row in gdf.iterrows():
        # ── derive a clean ROI name ──────────────────────────────────────
        props = row.get("properties") or {}
        if isinstance(props, str):          # geopandas sometimes flattens
            try:
                props = json.loads(props)
            except Exception:
                props = {}

        name = (
            props.get("name")
            or props.get("Name")
            or (props.get("classification") or {}).get("name")
            or props.get("objectType")
            or f"ROI_{i + 1}"
        )
        name = str(name).strip().replace(" ", "_")

        # ── spatial filter ───────────────────────────────────────────────
        mask = points.within(row.geometry)
        n_cells = mask.sum()
        print(f"[step 2]   {name}: {n_cells} cells")

        if n_cells == 0:
            print(f"[warn]     No cells found in {name} — skipping.")
            continue

        sub = adata[mask.values].copy()
        sub.uns["roi_name"] = name
        roi_adatas[name] = sub

    if not roi_adatas:
        raise RuntimeError(
            "No cells were assigned to any ROI. "
            "Check that the GeoJSON coordinates match the Xenium pixel space."
        )

    return roi_adatas

# ════════════════════════════════════════════════════════════════════════════
# 3 ─ Write one H5AD per ROI
# ════════════════════════════════════════════════════════════════════════════

def write_h5ads(roi_adatas: dict[str, ad.AnnData],
                output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    for name, sub in roi_adatas.items():
        out_path = output_dir / f"{name}.h5ad"
        sub.write_h5ad(str(out_path))
        print(f"[step 3] Written: {out_path}  ({sub.n_obs} cells)")
        written[name] = out_path

    return written


# ════════════════════════════════════════════════════════════════════════════
# 4 ─ Write metadata TSV
# ════════════════════════════════════════════════════════════════════════════

def write_tsv(roi_adatas: dict[str, ad.AnnData],
              written_paths: dict[str, Path],
              data_dir: Path,
              output_dir: Path) -> None:
    """
    Columns: sample_id | fastq_dir | slide | area

    fastq_dir  →  data_dir itself (Xenium outputs contain FASTQ-equivalent
                  transcriptomics data; adjust if you have a separate FASTQ
                  folder).
    """
    tsv_path = output_dir / "metadata.tsv"
    fastq_dir = str(data_dir)           # adjust if your FASTQs live elsewhere

    lines = ["sample_id\tfastq_dir\tslide\tarea"]
    for name in roi_adatas:
        lines.append(f"{name}\t{fastq_dir}\t{SLIDE_ID}\t{AREA_ID}")

    tsv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[step 4] TSV written: {tsv_path}")


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def main() -> None:
    # 0 ─ ensure data exists
    ensure_data(DATA_DIR)

    # 1 ─ build AnnData
    adata = load_xenium_anndata(DATA_DIR)

    # 2 ─ find GeoJSON and split
    geojson_path = find_geojson(DATA_DIR)
    print(f"[step 2] Using GeoJSON: {geojson_path.name}")
    roi_adatas = split_by_roi(adata, geojson_path, data_dir=DATA_DIR)

    # 3 ─ write H5ADs
    written_paths = write_h5ads(roi_adatas, OUTPUT_DIR)

    # 4 ─ write TSV
    write_tsv(roi_adatas, written_paths, DATA_DIR, OUTPUT_DIR)

    print("\n✓ All done.")
    print(f"  H5ADs  → {OUTPUT_DIR}")
    print(f"  TSV    → {OUTPUT_DIR / 'metadata.tsv'}")


if __name__ == "__main__":
    main()
