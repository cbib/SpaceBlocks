#!/usr/bin/env python3
"""
extract_external_annotation.py
==============================
Extract cell-type annotations from already-annotated .h5ad files into per-sample
metadata TSVs that SpaceBlocks' external-annotation feature can consume.

For each input .h5ad it writes ``{outdir}/metadata_{sample}.tsv``: a tab-separated
table indexed by cell barcode (``adata.obs_names``) containing the requested
annotation column(s). To then test the feature, point
``config["precomputed_metadata_dir"]`` at ``{outdir}`` and set::

    external_annotation:
      enabled: true
      column: "<out-column>"   # the column name written below

`annotate_cells` reads that column from the metadata (preferring
precomputed_metadata_dir over the pipeline's own metadata_{sample}.tsv), reindexes it
to the current cells (unmatched → "Unannotated"), and writes obs["cell_type_external"].

NOTE: the barcodes (obs_names) in the annotated files must match the pipeline's
adata.obs_names for the reindex to line up; otherwise cells fall back to "Unannotated".

Usage
-----
    python tools/extract_external_annotation.py \
        /scratch/.../adata_S1_annotated.h5ad /scratch/.../adata_S2_annotated.h5ad \
        --column manual_celltype --outdir .tests/external_metadata

    # sample name from an obs column instead of the filename:
    python tools/extract_external_annotation.py /scratch/*.h5ad \
        --column manual_celltype --sample-obs sample --outdir OUT
"""

import argparse
import os
import re
import sys

import anndata as ad


def derive_sample(path, regex):
    """Sample id from the filename: a --sample-regex capture group, else the file
    stem with the common `adata_` prefix / `_annotated` suffix stripped."""
    stem = os.path.splitext(os.path.basename(path))[0]
    if regex:
        m = re.search(regex, stem)
        if not m:
            sys.exit(f"[error] --sample-regex {regex!r} did not match '{stem}'")
        return m.group(1)
    stem = re.sub(r"^adata_", "", stem)
    stem = re.sub(r"_annotated$", "", stem)
    return stem


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5ads", nargs="+", help="annotated .h5ad file(s)")
    ap.add_argument("--column", "-c", action="append", required=True,
                    help="obs column holding the annotation (repeatable)")
    ap.add_argument("--outdir", "-o", required=True, help="output directory")
    ap.add_argument("--out-column", default=None,
                    help="rename the (single) column on write; default keeps the name")
    ap.add_argument("--sample-regex", default=None,
                    help="regex with ONE capture group to derive the sample from the filename")
    ap.add_argument("--sample-obs", default=None,
                    help="obs column holding the sample name (overrides filename derivation)")
    args = ap.parse_args()

    if args.out_column and len(args.column) != 1:
        sys.exit("[error] --out-column requires exactly one --column")
    os.makedirs(args.outdir, exist_ok=True)

    for path in args.h5ads:
        if not os.path.isfile(path):
            sys.exit(f"[error] no such file: {path}")
        adata = ad.read_h5ad(path, backed="r")   # obs only; X is not loaded
        obs = adata.obs
        missing = [c for c in args.column if c not in obs.columns]
        if missing:
            sys.exit(f"[error] column(s) {missing} not in obs of {path}.\n"
                     f"        available obs columns: {list(obs.columns)}")

        if args.sample_obs:
            if args.sample_obs not in obs.columns:
                sys.exit(f"[error] --sample-obs '{args.sample_obs}' not in obs of {path}")
            vals = list(dict.fromkeys(obs[args.sample_obs].astype(str)))
            if len(vals) != 1:
                sys.exit(f"[error] {path}: obs['{args.sample_obs}'] has {len(vals)} distinct "
                         "values; expected one sample per file")
            sample = vals[0]
        else:
            sample = derive_sample(path, args.sample_regex)

        out = obs[list(args.column)].copy()
        out.index.name = "barcode"
        if args.out_column:
            out.columns = [args.out_column]

        dst = os.path.join(args.outdir, f"metadata_{sample}.tsv")
        out.to_csv(dst, sep="\t")
        try:
            adata.file.close()
        except Exception:
            pass
        print(f"{os.path.basename(path)} -> {dst}  "
              f"({len(out)} cells; {out.columns[0]}: {out.iloc[:, 0].nunique()} types)")


if __name__ == "__main__":
    main()
