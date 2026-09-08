#!/usr/bin/env python3
"""
normalize_external_metadata.py — normalise SpaceBlocks external-annotation metadata.

Rewrites metadata_{sample}.tsv into the strict 2-column shape annotate_cells expects:

    barcode <TAB> <label-column>
    <plain cell_id matching adata.obs_names> <TAB> <label>

Ragged files (unnamed suffixed index + a real cell_id column + label) are re-keyed on
their plain cell_id column; already-correct 2-column files pass through unchanged
(the function is idempotent).

Two entry points:
  * normalize_sample_metadata(...) — called by the head blocks (prepare_input_x5k /
    prepare_input_ate) for their one sample, with the directory taken from
    config["precomputed_metadata_dir"]. Runs per sample, so parallel head jobs never
    touch the same file.
  * a small __main__ CLI, for a one-off pass over a whole directory.

Files are rewritten IN PLACE in the configured metadata directory. The untouched
original is preserved (write-once) in a sibling folder named 'original_metadata', so
the very first raw copy is always recoverable even across re-runs.
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys

import pandas as pd


def original_backup_dir(metadata_dir: str) -> str:
    """Sibling 'original_metadata' folder next to the configured metadata dir."""
    md = metadata_dir.rstrip("/")
    return os.path.join(os.path.dirname(md), "original_metadata")


def _normalise_frame(df: pd.DataFrame, label_col: str, id_col: str, src: str) -> pd.DataFrame:
    """Return a 2-column frame indexed by the plain cell_id ('barcode') + the label."""
    if label_col not in df.columns:
        raise ValueError(
            f"label column '{label_col}' not found in {src}; columns seen: {list(df.columns)}"
        )
    if id_col in df.columns:
        key = df[id_col].astype(str).values        # ragged file: real cell_id in a column
    else:
        key = df.index.astype(str).values          # already-correct: index is the cell_id
    out = pd.DataFrame({label_col: df[label_col].astype(str).values})
    out.index = key
    out.index.name = "barcode"
    if out.index.duplicated().any():
        n = int(out.index.duplicated().sum())
        print(f"  [warn] {os.path.basename(src)}: {n} duplicate barcode(s); keeping first")
        out = out[~out.index.duplicated(keep="first")]
    return out


def normalize_sample_metadata(
    metadata_dir: str,
    sample_id: str,
    label_column: str = "celltype_annotation",
    id_column: str = "cell_id",
    logger=None,
):
    """
    Normalise metadata_{sample_id}.tsv inside `metadata_dir`, in place.

    Before the first rewrite the untouched original is copied to the sibling
    'original_metadata' folder (write-once: later runs never clobber it). Reading the
    current live file every time keeps this idempotent — a file that is already in the
    2-column shape passes straight through.

    Returns the path rewritten, or None when there is no metadata file for this sample.
    """
    log = (logger.info if logger is not None else print)
    warn = (logger.warning if logger is not None else print)

    if not metadata_dir:
        return None
    fname = f"metadata_{sample_id}.tsv"
    live = os.path.join(metadata_dir, fname)
    if not os.path.isfile(live):
        log(f"[normalize] no external metadata for '{sample_id}' ({live}); skipping.")
        return None

    backup_dir = original_backup_dir(metadata_dir)
    os.makedirs(backup_dir, exist_ok=True)          # concurrent-safe across parallel jobs
    backup = os.path.join(backup_dir, fname)
    if not os.path.exists(backup):
        shutil.copy2(live, backup)                  # preserve the first-seen original
        log(f"[normalize] backed up original → {backup}")

    df = pd.read_csv(live, sep="\t", index_col=0, comment="#", dtype=str)
    out = _normalise_frame(df, label_column, id_column, live)
    out.to_csv(live, sep="\t")                      # rewrite in the configured dir
    log(f"[normalize] {fname}: {len(out)} cells, {out[label_column].nunique()} types "
        f"→ {live}")
    return live


def _cli():
    ap = argparse.ArgumentParser(description="Normalise external-annotation metadata TSVs.")
    ap.add_argument("--metadata-dir", required=True,
                    help="directory of metadata_*.tsv (== config precomputed_metadata_dir)")
    ap.add_argument("--label-column", default="celltype_annotation")
    ap.add_argument("--id-column", default="cell_id")
    ap.add_argument("--sample", default=None,
                    help="normalise only metadata_<sample>.tsv (default: every file in dir)")
    args = ap.parse_args()

    if args.sample:
        normalize_sample_metadata(args.metadata_dir, args.sample,
                                  args.label_column, args.id_column)
        return

    files = sorted(glob.glob(os.path.join(args.metadata_dir, "metadata_*.tsv")))
    if not files:
        sys.exit(f"[error] no metadata_*.tsv in {args.metadata_dir}")
    print(f"Normalising {len(files)} file(s) in {args.metadata_dir}:")
    for f in files:
        s = os.path.basename(f)[len("metadata_"):-len(".tsv")]
        normalize_sample_metadata(args.metadata_dir, s, args.label_column, args.id_column)


if __name__ == "__main__":
    _cli()
