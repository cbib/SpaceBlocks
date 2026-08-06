"""
contract_io.py — shared dtype hygiene for the contract h5ad.
============================================================
Newer pandas/anndata can hand back pandas EXTENSION dtypes (StringDtype indices,
nullable Boolean/Int64/Float64 columns) where scanpy/scipy expect plain numpy. The
classic symptom is ``adata.var_names.str.startswith(...)`` returning a nullable
BooleanArray, which then makes ``sc.pp.calculate_qc_metrics`` raise
``"'BooleanArray' object has no attribute 'nonzero'"``.

This lived inline in a single head (prepare_input_ate) before, which left decoupled
contracts and the other core readers unprotected. It is now one helper, applied at
BOTH boundaries so every input is covered exactly once, idempotently:

  * the head writes a clean contract (prepare_input_*), and
  * every core reader of the RAW contract (preprocess_umap, qc_sweep) normalises on
    load — which also covers decoupled contracts supplied directly, and any future
    head that forgets to.

``validate_input`` only inspects the contract (and tolerates extension dtypes), so it
deliberately does not call this.
"""
import logging

import pandas as pd

_log = logging.getLogger("contract_io")


def normalize_contract_dtypes(adata, logger=None):
    """Coerce pandas EXTENSION dtypes in obs/var (and the axis indices) to plain numpy.

    Categoricals are preserved. Mutates ``adata`` in place and returns it. Idempotent:
    a second call finds nothing left to convert.
    """
    log = logger or _log

    def _plain_index(idx):
        return pd.Index(idx.astype(object), dtype=object)

    adata.obs_names = _plain_index(adata.obs_names)
    adata.var_names = _plain_index(adata.var_names)

    converted = []
    for frame_name, frame in (("obs", adata.obs), ("var", adata.var)):
        for col in list(frame.columns):
            dt = frame[col].dtype
            if isinstance(dt, pd.CategoricalDtype):
                continue                       # categoricals are wanted, keep them
            if not pd.api.types.is_extension_array_dtype(dt):
                continue
            n_na = int(frame[col].isna().sum())
            if isinstance(dt, pd.StringDtype):
                frame[col] = frame[col].astype(object)
            elif isinstance(dt, pd.BooleanDtype):
                if n_na:
                    log.warning("  %s['%s']: %d missing value(s) in a nullable boolean "
                                "column set to False during conversion.",
                                frame_name, col, n_na)
                frame[col] = frame[col].fillna(False).astype(bool)
            else:                              # nullable Int64/Float64 and friends
                frame[col] = frame[col].astype(
                    "float64" if n_na else frame[col].dtype.numpy_dtype)
            converted.append(f"{frame_name}['{col}'] {dt}")
    if converted:
        log.info("Normalised %d pandas extension column(s) to numpy dtypes: %s",
                 len(converted), ", ".join(converted))
    return adata
