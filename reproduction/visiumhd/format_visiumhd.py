#!/usr/bin/env python3
"""Visium HD public dataset -> 3 synthetic SpaceBlocks samples (mode: decoupled).

Strategy (deliberately cheap, but analysis-ready):
  1. load one bin resolution and keep the TOP_GENES most-expressed genes;
  2. assign every bin a region from the QuPath GeoJSON polygons
     (obs["region_annotation"] = "Region 1".."Region N"; bins outside any -> "Unlabeled");
  3. randomly split all bins into N_SAMPLES non-overlapping samples, balanced *within* each
     region (seeded). Every sample then carries every region, so integrate_samples has several
     samples to align and pseudobulk_de has region replicates for DESeq2.

Each sample is written as a contract h5ad (X = raw counts, obsm["spatial"],
obs = sample / cell_id / region_annotation), plus a core_samples.tsv.

Download the dataset manually (see reproduction.md) into ./data/<SAMPLE>/, then run this.
deps: numpy pandas h5py scipy anndata geopandas shapely pyarrow
"""
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
import anndata as ad
import geopandas as gpd
from shapely.geometry import Point

SAMPLE    = "Visium_HD_Mouse_Brain"        # 10x dataset prefix (must match your GeoJSON slide)
BIN_SIZE  = "008um"                        # 002um | 008um | 016um
DATA_DIR  = Path("data") / SAMPLE          # the 10x files live here (manual download)
GEOJSON   = Path("demo_vhd.geojson")       # QuPath export on the hires image
OUT_DIR   = Path("contracts")              # -> contracts/<sampleK>/<sampleK>_unfiltered.h5ad
TOP_GENES = 500                            # keep only the N most-expressed genes (light test)
N_SAMPLES = 3                              # random non-overlapping samples to split into
SEED      = 0                              # reproducible split


def load_bins(data_dir: Path) -> ad.AnnData:
    """Read one bin resolution into AnnData with full-res pixel coordinates."""
    for tar in (f"{SAMPLE}_binned_outputs.tar.gz", f"{SAMPLE}_spatial.tar.gz"):
        sub = data_dir / tar.replace(f"{SAMPLE}_", "").replace(".tar.gz", "")
        if (data_dir / tar).exists() and not sub.exists():
            with tarfile.open(data_dir / tar) as tf:
                tf.extractall(data_dir)
    bin_dir = data_dir / "binned_outputs" / f"square_{BIN_SIZE}"
    with h5py.File(bin_dir / "filtered_feature_bc_matrix.h5") as f:
        m = f["matrix"]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                          shape=m["shape"][:]).T.tocsr()
        barcodes = m["barcodes"][:].astype(str)
        var = pd.DataFrame(index=m["features"]["name"][:].astype(str))
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=var)
    a.var_names_make_unique()
    pos = pd.read_parquet(bin_dir / "spatial" / "tissue_positions.parquet")
    pos = pos.set_index("barcode").loc[a.obs_names]
    a.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
    a.uns["hires_scalef"] = json.load(
        open(bin_dir / "spatial" / "scalefactors_json.json"))["tissue_hires_scalef"]
    return a


def top_genes(a: ad.AnnData, n: int) -> ad.AnnData:
    """Keep the n most-expressed genes (shared across all bins)."""
    totals = np.asarray(a.X.sum(axis=0)).ravel()
    return a[:, np.sort(np.argsort(totals)[::-1][:n])].copy()


def assign_regions(a: ad.AnnData, geojson: Path) -> None:
    """obs['region_annotation'] = 'Region i' per polygon (QuPath annotates the hires image)."""
    gdf = gpd.read_file(geojson)
    gdf["geometry"] = gdf.scale(xfact=1.0 / a.uns["hires_scalef"],
                                yfact=1.0 / a.uns["hires_scalef"], origin=(0, 0))
    pts = gpd.GeoSeries([Point(*xy) for xy in a.obsm["spatial"]])
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
    a = load_bins(DATA_DIR)
    a = top_genes(a, TOP_GENES)
    assign_regions(a, GEOJSON)
    samples = make_samples(a, N_SAMPLES, SEED)
    for name, sub in samples.items():
        out = OUT_DIR / f"{name}.h5ad"   # flat: decoupled reads <dir>/<sample>.h5ad
        out.parent.mkdir(parents=True, exist_ok=True)
        sub.write_h5ad(out)
        print(f"  {name}: {sub.n_obs} bins")
    (OUT_DIR / "core_samples.tsv").write_text("sample\n" + "\n".join(samples) + "\n")
    print(f"[done] {len(samples)} samples x {a.n_vars} genes -> {OUT_DIR}/  (+ core_samples.tsv)")


if __name__ == "__main__":
    main()
