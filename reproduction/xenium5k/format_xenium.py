#!/usr/bin/env python3
"""Xenium 5K public dataset -> 3 synthetic SpaceBlocks samples (mode: decoupled).

Strategy (deliberately cheap, but analysis-ready):
     demo marker sets in gene_queries_demo.tsv survive the trim);
  2. embed a greyscale morphology composite in uns["spatial"], so the CoreBlock's
     spatial plots have a tissue background instead of a bare scatter;
  3. assign every cell a region from the QuPath GeoJSON polygons
     (obs["region_annotation"] = "Region 1".."Region N"; cells outside any -> "Unlabeled");
  4. randomly split all cells into N_SAMPLES non-overlapping samples, balanced *within* each
     region (seeded). Every sample then carries every region, so integrate_samples has several
     samples to align and pseudobulk_de has region replicates for DESeq2 dispersion estimation.

Each sample is written as a contract h5ad (X = raw counts, obsm["spatial"], uns["spatial"],
obs = sample / cell_id / region_annotation), plus a core_samples.tsv.

Download the bundle manually (see reproduction.md) into ./data/<SAMPLE>/, then run this.
deps: numpy pandas anndata geopandas shapely spatialdata-io  (envs/xenium5k.yaml)
"""
import json
from pathlib import Path

import numpy as np
import anndata as ad
import geopandas as gpd
from shapely.geometry import Point
from shapely.affinity import scale as shp_scale

SAMPLE     = "Xenium_Prime_Human_Skin_FFPE"
DATA_DIR   = Path("data") / SAMPLE         # Xenium output bundle (manual download)
GEOJSON    = Path("demo_x5k.geojson")      # QuPath export on the morphology image
OUT_DIR    = Path("contracts")             # -> contracts/<sampleK>.h5ad (flat, decoupled)
TOP_GENES  = 500                           # keep only the N most-expressed genes (light test)
N_SAMPLES  = 3                             # random non-overlapping samples to split into
SEED       = 0                             # reproducible split
HIRES_LEVEL = 3                            # morphology pyramid level embedded as the background


def load_xenium(data_dir: Path):
    """Return (AnnData with micron centroids, SpatialData object)."""
    import spatialdata_io
    sdata = spatialdata_io.xenium(str(data_dir), cells_as_circles=True)
    a = sdata.tables[list(sdata.tables)[0]]
    if "spatial" not in a.obsm:
        c = sdata["cell_circles"]
        a.obsm["spatial"] = np.column_stack([c.geometry.x, c.geometry.y]).astype(np.float32)
    return a, sdata


def top_genes(a: ad.AnnData, n: int) -> ad.AnnData:
    """Keep the n most-expressed genes (shared across all cells)."""
    totals = np.asarray(a.X.sum(axis=0)).ravel()
    return a[:, np.sort(np.argsort(totals)[::-1][:n])].copy()


def _pixel_size(data_dir: Path, default: float = 0.2125) -> float:
    for p in (data_dir / "experiment.xenium", data_dir / "experiment.xenium.json"):
        if p.exists():
            meta = json.load(open(p))
            for k in ("pixel_size", "pixelSize", "micronsPerPixel"):
                if k in meta:
                    return float(meta[k])
    return default


def grey_composite(sdata, level: int):
    """Greyscale mean-projection of the morphology_focus channels at a pyramid level.

    Mirrors workflow/scripts/xenium_image.py (used by the real Xenium head): channel-agnostic,
    so it works whatever the channel identities/order are. Returns (uint8 RGB, downsample factor).
    """
    if "morphology_focus" not in sdata.images:
        return None, None
    img = sdata.images["morphology_focus"]
    scales = list(img.keys())
    ds = img[scales[min(level, len(scales) - 1)]]
    arr = ds[list(ds.data_vars)[0]].values
    if hasattr(arr, "compute"):
        arr = arr.compute()
    ds0 = img[scales[0]]
    full_h = ds0[list(ds0.data_vars)[0]].shape[1]
    norm = []
    for ch in arr:                                  # per-channel percentile normalisation
        lo, hi = np.percentile(ch, (1, 99))
        norm.append(np.clip((ch - lo) / max(hi - lo, 1e-9), 0, 1))
    grey = np.mean(norm, axis=0)
    rgb = (np.stack([grey, grey, grey], axis=-1) * 255).astype(np.uint8)
    return rgb, arr.shape[1] / full_h


def assign_regions(a: ad.AnnData, geojson: Path, px: float) -> None:
    """obs['region_annotation'] = 'Region i' per polygon (cells outside any -> Unlabeled)."""
    gdf = gpd.read_file(geojson)
    xy = a.obsm["spatial"]
    # QuPath annotates in pixels; Xenium cells are in microns. Rescale if footprints differ.
    gspan = gdf.total_bounds[2] - gdf.total_bounds[0]
    cspan = xy[:, 0].max() - xy[:, 0].min()
    if cspan and not (1 / 1.5 < gspan / cspan < 1.5):
        gdf["geometry"] = gdf.geometry.apply(
            lambda g: shp_scale(g, xfact=px, yfact=px, origin=(0, 0)))
    pts = gpd.GeoSeries([Point(*p) for p in xy])
    region = np.array(["Unlabeled"] * a.n_obs, dtype=object)
    for i, row in gdf.iterrows():
        region[pts.within(row.geometry).values] = f"Region {i + 1}"
    a.obs["region_annotation"] = region


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
    a, sdata = load_xenium(DATA_DIR)
    a = top_genes(a, TOP_GENES)
    px = _pixel_size(DATA_DIR)
    rgb, ds_factor = grey_composite(sdata, HIRES_LEVEL)
    if rgb is None:
        print("[warn] no morphology_focus image — spatial plots will be a bare scatter")
    assign_regions(a, GEOJSON, px)
    samples = make_samples(a, N_SAMPLES, SEED)
    for name, sub in samples.items():
        if rgb is not None:
            # micron -> level pixel:  micron / pixel_size_um * downsample_factor
            sub.uns["spatial"] = {name: {
                "images": {"hires": rgb},
                "scalefactors": {"tissue_hires_scalef": float(ds_factor / px),
                                 "spot_diameter_fullres": 10.0,
                                 "pixel_size_um": px}}}
        out = OUT_DIR / f"{name}.h5ad"   # flat: decoupled reads <dir>/<sample>.h5ad
        out.parent.mkdir(parents=True, exist_ok=True)
        sub.write_h5ad(out)
        print(f"  {name}: {sub.n_obs} cells")
    (OUT_DIR / "core_samples.tsv").write_text("sample\n" + "\n".join(samples) + "\n")
    img_note = f", image {rgb.shape[1]}x{rgb.shape[0]} px" if rgb is not None else ""
    print(f"[done] {len(samples)} samples x {a.n_vars} genes{img_note} -> {OUT_DIR}/")


if __name__ == "__main__":
    main()
