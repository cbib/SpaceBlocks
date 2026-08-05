#!/usr/bin/env python3
"""Regenerate the deterministic synthetic contract h5ads used by CI and catalog checks."""

from pathlib import Path

import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

SEED = 0
SAMPLES = ("S1", "S2")
N_CELLS = 60
GENES = [f"GENE{i:02d}" for i in range(28)] + ["MT-CO1", "MT-ND1"]
REGIONS = ("Tumor area", "Healthy area")
TEST_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = TEST_ROOT / "data"


def build_sample(sample: str, rng: np.random.Generator) -> ad.AnnData:
    """Build one tiny AnnData object satisfying the SpaceBlocks input contract."""
    matrix = csr_matrix(
        rng.poisson(1.0, size=(N_CELLS, len(GENES))).astype("int64")
    )
    cell_ids = [f"{sample}_{i}" for i in range(N_CELLS)]
    adata = ad.AnnData(
        X=matrix,
        obs={
            "sample": [sample] * N_CELLS,
            "cell_id": cell_ids,
            "region_annotation": rng.choice(REGIONS, size=N_CELLS),
        },
    )
    adata.obs_names = cell_ids
    adata.var_names = GENES
    adata.obsm["spatial"] = rng.uniform(0, 1000, size=(N_CELLS, 2))
    return adata


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    for sample in SAMPLES:
        output = DATA_DIR / f"{sample}.h5ad"
        adata = build_sample(sample, rng)
        adata.write_h5ad(output)
        print(f"wrote {output} {adata.shape}")


if __name__ == "__main__":
    main()
