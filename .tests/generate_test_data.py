#!/usr/bin/env python3
"""Generate tiny synthetic *contract* h5ads for the CI smoke test (mode: decoupled).

Each file matches the contract: raw integer counts in X, obsm['spatial'], and obs
columns sample / cell_id / region_annotation. Deliberately minimal (~60 cells, 30
genes) so the repo stays catalog-clonable and CI is fast.
"""
import os
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix

RNG = np.random.default_rng(0)
GENES = [f"GENE{i:02d}" for i in range(28)] + ["MT-CO1", "MT-ND1"]  # 2 mito genes
REGIONS = ["Tumor area", "Healthy area"]
HERE = os.path.dirname(os.path.abspath(__file__))

for sample in ("S1", "S2"):
    n = 60
    X = csr_matrix(RNG.poisson(1.0, size=(n, len(GENES))).astype("int64"))
    obs = dict(
        sample=[sample] * n,
        cell_id=[f"{sample}_{i}" for i in range(n)],
        region_annotation=RNG.choice(REGIONS, size=n),
    )
    a = ad.AnnData(X=X, obs=obs)
    a.obs_names = [f"{sample}_{i}" for i in range(n)]
    a.var_names = GENES
    a.obsm["spatial"] = RNG.uniform(0, 1000, size=(n, 2))
    out = os.path.join(HERE, "integration", "data", sample, f"{sample}_unfiltered.h5ad")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    a.write_h5ad(out)
    print("wrote", out, a.shape)
