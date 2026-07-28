#!/usr/bin/env python3
"""Visium HD (Space Ranger >= 4.0) public dataset -> 3 synthetic SpaceBlocks samples.

Space Ranger >= 4.0 adds H&E nucleus/cell segmentation, so Visium HD yields real single cells
(not 2-16 um bins). This reads the segmented outputs **the same way the Visium HD HeadBlock does**
(scanpy + a one-line image read), so no SpatialData/zarr intermediate is needed:

  1. load the segmented cell matrix (raw counts) + per-cell full-res pixel centroids
     (mean of each cell's 2 um bins, from barcode_mappings + tissue_positions — head-identical);
  2. embed the hires H&E image + native scalefactor in uns["spatial"] (tissue background);
  3. assign every cell a region from the QuPath GeoJSON polygons
     (obs["region_annotation"] = "Region 1".."Region N"; cells outside any -> "Unlabeled");
  4. keep the TOP_GENES most-expressed genes, then randomly split all cells into N_SAMPLES
     non-overlapping samples, balanced *within* each region (seeded) — so integrate_samples has
     several samples and pseudobulk_de has region replicates for DESeq2.

Each sample is a flat contract h5ad (X = raw counts, obsm["spatial"], uns["spatial"] with the
hires H&E, obs = sample / cell_id / region_annotation), plus a core_samples.tsv. Coordinates are
stored in HIRES-image pixel space (embedded scalefactor = 1.0): the CoreBlock's spatial plots use
a fixed spot_size, which renders invisibly against full-res coordinates, so the demo remaps to the
image's own pixel space. The Visium HD HeadBlock keeps full-res coordinates + the native
scalefactor; this remapping is local to the demo contract.

Needs (all from the 10x page): segmented_outputs/, barcode_mappings.parquet, and the
binned_outputs spatial/tissue_positions.parquet. See reproduction.md.
deps: numpy pandas anndata scanpy pillow geopandas shapely  (envs/visiumhd.yaml)
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import geopandas as gpd
from PIL import Image
from shapely.geometry import Point
from shapely.affinity import scale as shp_scale

SAMPLE      = "Visium_HD_Mouse_Brain"      # 10x dataset prefix (must match your GeoJSON slide)
DATA_DIR    = Path("data") / SAMPLE        # the Space Ranger outs live here (manual download)
GEOJSON     = Path("demo_vhd.geojson")      # QuPath export on the hires image
OUT_DIR     = Path("contracts")            # -> contracts/<sampleK>.h5ad (flat, decoupled)
TOP_GENES   = 500                          # keep only the N most-expressed genes (light test)
N_SAMPLES   = 3                            # random non-overlapping samples to split into
SEED        = 0                            # reproducible split


def _find(*relpaths):
    """First existing path under DATA_DIR among the given candidates (SR layout varies)."""
    for rel in relpaths:
        p = DATA_DIR / rel
        if p.exists():
            return p
    raise FileNotFoundError(f"none of these exist under {DATA_DIR}: {relpaths}")


def _norm_id(raw):
    """Canonical integer-string cell id, matching the Visium HD head's _normalise_cell_id:
    'cellid_000000123-1' -> '123', 'cellid_...' -> digits, plain numeric kept."""
    s = str(raw)
    if s.endswith("-1"):
        s = s[:-2]
    if s.lower().startswith("cellid_"):
        s = s.split("_", 1)[1]
    s = s.lstrip("0") or "0"
    return s


def load_cells() -> ad.AnnData:
    """Segmented cell matrix + full-res pixel centroids + hires H&E image.

    Mirrors the Visium HD HeadBlock exactly: cell centroids are the mean full-res pixel position
    of each cell's constituent 2 um bins, from barcode_mappings.parquet + tissue_positions.parquet.
    This is unambiguous (no GeoJSON/CRS guessing) and puts obsm["spatial"] in the same full-res
    pixel space the native tissue_hires_scalef expects. Needs, in addition to segmented_outputs/,
    the top-level barcode_mappings.parquet and the binned_outputs spatial/tissue_positions.parquet.
    """
    seg = _find("outs/segmented_outputs", "segmented_outputs")
    mtx = _find(str(seg.relative_to(DATA_DIR) / "filtered_feature_cell_matrix.h5"),
                str(seg.relative_to(DATA_DIR) / "raw_feature_cell_matrix.h5"))
    a = sc.read_10x_h5(str(mtx))
    a.var_names_make_unique()
    a.obs["_cid"] = a.obs_names.map(_norm_id)

    # bin -> cell mapping and bin full-res pixel positions
    bmap = pd.read_parquet(_find("Visium_HD_Mouse_Brain_barcode_mappings.parquet",
                                 "barcode_mappings.parquet", "outs/barcode_mappings.parquet"))
    tpos = pd.read_parquet(_find(
        "outs/binned_outputs/square_002um/spatial/tissue_positions.parquet",
        "binned_outputs/square_002um/spatial/tissue_positions.parquet",
        "outs/segmented_outputs/spatial/tissue_positions.parquet",
        "segmented_outputs/spatial/tissue_positions.parquet")).set_index("barcode")

    cid_col = next(c for c in bmap.columns if c.lower() == "cell_id")
    bin_col = next(c for c in bmap.columns if c.lower() == "square_002um")
    bmap["_cid"] = bmap[cid_col].astype(str).map(_norm_id)
    bmap = bmap[bmap["_cid"].isin(a.obs["_cid"]) & (bmap["_cid"] != "0")]
    bmap = bmap.merge(tpos[["pxl_row_in_fullres", "pxl_col_in_fullres"]],
                      left_on=bin_col, right_index=True, how="inner")
    cent = bmap.groupby("_cid")[["pxl_col_in_fullres", "pxl_row_in_fullres"]].mean()

    has = a.obs["_cid"].isin(cent.index)
    if has.sum() == 0:
        raise RuntimeError(
            f"no cell matched between matrix and barcode_mappings. matrix ids "
            f"e.g. {list(a.obs['_cid'][:3])}; centroid ids e.g. {list(cent.index[:3])}")
    a = a[has].copy()
    fullres_xy = cent.loc[a.obs["_cid"].values].to_numpy(np.float64)          # (x=col, y=row)
    del a.obs["_cid"]

    # hires H&E image + native scalefactor + real spot diameter — the head's one-liner.
    img_path = _find(str(seg.relative_to(DATA_DIR) / "spatial/tissue_hires_image.png"),
                     "outs/spatial/tissue_hires_image.png", "spatial/tissue_hires_image.png")
    sf = json.load(open(_find(str(seg.relative_to(DATA_DIR) / "spatial/scalefactors_json.json"),
                              "outs/spatial/scalefactors_json.json",
                              "spatial/scalefactors_json.json")))
    hires_scalef = sf["tissue_hires_scalef"]
    # Demo convention: store coordinates in HIRES-image pixel space (full-res * hires_scalef) and
    # embed the image with scalef 1.0. scanpy sizes spots by spot_size * scalef, so full-res coords
    # (which force scalef ~0.13) make the CoreBlock's hardcoded spot_size=20 render invisibly; in
    # hires space (scalef 1.0) the same spot_size=20 is clearly visible. The Visium HD HeadBlock
    # keeps full-res coordinates + the native scalefactor; this remapping is local to the demo.
    a.obsm["spatial"] = fullres_xy * hires_scalef
    a.uns["_hires_img"] = np.asarray(Image.open(img_path))
    a.uns["_hires_scalef_orig"] = hires_scalef      # kept only to reconcile the region GeoJSON
    return a


def top_genes(a: ad.AnnData, n: int) -> ad.AnnData:
    totals = np.asarray(a.X.sum(axis=0)).ravel()
    keep = a[:, np.sort(np.argsort(totals)[::-1][:n])].copy()
    keep.uns["_hires_img"] = a.uns["_hires_img"]        # carry through the gene subset
    keep.uns["_hires_scalef_orig"] = a.uns["_hires_scalef_orig"]
    return keep


def assign_regions(a: ad.AnnData, geojson: Path, hires_scalef: float) -> None:
    """obs['region_annotation'] = 'Region i' per polygon (cells outside any -> Unlabeled).

    Cells live in FULL-RES pixels. QuPath annotations are usually drawn on the hires image
    (tissue_hires_image.png), so their polygons are in HIRES pixels and must be scaled up by
    1/hires_scalef to reach full-res. We auto-detect which space the GeoJSON is in by comparing
    its span to the cells' span, snapping to the known factor {1, 1/hires_scalef} — so it works
    whether you annotated the hires image or the full-res image."""
    gdf = gpd.read_file(geojson)
    xy = a.obsm["spatial"]
    gspan = max(gdf.total_bounds[2] - gdf.total_bounds[0], gdf.total_bounds[3] - gdf.total_bounds[1])
    cspan = max(xy[:, 0].max() - xy[:, 0].min(), xy[:, 1].max() - xy[:, 1].min())
    ratio = cspan / gspan if gspan else 1.0
    candidates = {hires_scalef: "coords>geojson", 1.0: "same space",
                  (1.0 / hires_scalef): "geojson>coords"}
    factor = min(candidates, key=lambda f: abs(f - ratio))
    if abs(factor - 1.0) > 1e-6:
        gdf["geometry"] = gdf.geometry.apply(
            lambda g: shp_scale(g, xfact=factor, yfact=factor, origin=(0, 0)))
    print(f"  region GeoJSON in {candidates[factor]} space "
          f"(span ratio {ratio:.2f} -> polygon scale x{factor:.3f})")
    pts = gpd.GeoSeries([Point(*p) for p in xy])
    region = np.array(["Unlabeled"] * a.n_obs, dtype=object)
    for i, row in gdf.iterrows():
        region[pts.within(row.geometry).values] = f"Region {i + 1}"
    a.obs["region_annotation"] = region
    vc = {k: int(v) for k, v in zip(*np.unique(region, return_counts=True))}
    print(f"  regions: {vc}")


def make_samples(a: ad.AnnData, n: int, seed: int) -> dict[str, ad.AnnData]:
    """Non-overlapping random split into n samples, balanced within each region."""
    rng = np.random.default_rng(seed)
    sid = np.empty(a.n_obs, dtype=object)
    reg = a.obs["region_annotation"].to_numpy()
    for r in np.unique(reg):
        idx = np.where(reg == r)[0]
        rng.shuffle(idx)
        for k, part in enumerate(np.array_split(idx, n)):
            sid[part] = f"sample{k + 1}"
    out = {}
    for s in sorted(set(sid)):
        sub = a[sid == s].copy()
        sub.obs["sample"] = s
        # contract convention: cell_id == obs_names (downstream TSVs are keyed by it)
        sub.obs["cell_id"] = sub.obs_names.astype(str)
        out[s] = sub
    return out


def main() -> None:
    a = load_cells()
    a = top_genes(a, TOP_GENES)
    img = a.uns.pop("_hires_img"); orig_scalef = a.uns.pop("_hires_scalef_orig")
    assign_regions(a, GEOJSON, orig_scalef)
    samples = make_samples(a, N_SAMPLES, SEED)
    for name, sub in samples.items():
        sub.uns["spatial"] = {name: {
            "images": {"hires": img},
            "scalefactors": {"tissue_hires_scalef": 1.0,   # coords are already in hires px
                             "spot_diameter_fullres": 1.0}}}
        out = OUT_DIR / f"{name}.h5ad"   # flat: decoupled reads <dir>/<sample>.h5ad
        out.parent.mkdir(parents=True, exist_ok=True)
        sub.write_h5ad(out)
        print(f"  {name}: {sub.n_obs} cells")
    (OUT_DIR / "core_samples.tsv").write_text("sample\n" + "\n".join(samples) + "\n")
    print(f"[done] {len(samples)} samples x {a.n_vars} genes, image {img.shape} -> {OUT_DIR}/")


if __name__ == "__main__":
    main()
