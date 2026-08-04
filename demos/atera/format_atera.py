#!/usr/bin/env python3
r"""
format_atera.py — turn the Atera supplemental CSVs into SpaceBlocks demo inputs.
===============================================================================
The 10x WTA preview ships two supplemental annotation files alongside the `outs/`
bundle. This script converts them into the formats SpaceBlocks consumes, so the demo
can exercise external annotation, neighbourhood analysis and spatial niches without
anyone hand-editing a 170,000-row table.

Nothing here is part of the workflow: it is a one-off preparation step, in the same
spirit as demos/visiumhd/format_visiumhd.py and demos/xenium5k/format_xenium.py.

Subcommands
-----------
metadata
    `<sample>_cell_groups.csv` (cell_id, group, color) -> `metadata_<sample>.tsv`,
    indexed by cell id with one label column. Point `precomputed_metadata_dir` at the
    output directory and set `external_annotation.enabled: true`. Also emits the
    matching `annotation_colors:` YAML block, using 10x's own display colours.

cluster-annotations
    Run AFTER `run_preprocessing`. Reads the pipeline's own metadata for the sample,
    labels each cell from the vendor annotation file, majority-votes those labels onto
    the Leiden clusters, and writes `cluster_annotations.tsv` in the layout
    `generate_annotation_template` produces (index `cluster`, one `<sample>_<res>`
    column). This resolves the ordering problem in the demo: `annotate_cells` needs a
    filled cluster_annotations TSV, but the clusters it refers to do not exist until
    preprocessing has run.

    The label source defaults to the external-annotation TSV written by `metadata`
    above, so the same file drives per-cell external annotation and the cluster
    mapping — one vocabulary, one place to correct it. Passing the raw
    `<sample>_cell_groups.csv` instead also works; the format is detected.

Usage
-----
    # 1. before running the pipeline
    python demos/atera/format_atera.py metadata \
        /scratch/.../WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv \
        --sample WTA_Preview_FFPE_Breast_Cancer \
        --outdir demos/atera/external_metadata

    # optional: confirm the cell ids line up with the contract before a long run
    python demos/atera/format_atera.py metadata ... --check-contract \
        /scratch/.../Samples/WTA_Preview_FFPE_Breast_Cancer/..._unfiltered.h5ad

    # 2. after `snakemake run_preprocessing`
    python demos/atera/format_atera.py cluster-annotations \
        --pipeline-metadata /scratch/.../Samples/<sample>/metadata_<sample>.tsv \
        --sample WTA_Preview_FFPE_Breast_Cancer --resolution 0.6
"""
import argparse
import os
import sys

import pandas as pd

LABEL_COLUMN = "vendor_celltype"


def read_cell_groups(path):
    """Read `<sample>_cell_groups.csv` -> DataFrame indexed by cell id.

    Cell ids are forced to str: Atera ids are alphanumeric, but a platform whose ids
    are all-digit would otherwise be inferred as int64 and silently fail to align with
    obs_names on the join."""
    df = pd.read_csv(path, dtype=str)
    missing = {"cell_id", "group"} - set(df.columns)
    if missing:
        sys.exit(f"[format_atera] {path} lacks column(s): {sorted(missing)}")
    df["cell_id"] = df["cell_id"].astype(str)
    df = df.set_index("cell_id")
    df["group"] = df["group"].fillna("Unannotated").str.strip()
    return df


def read_labels(path):
    """Read a per-cell label Series from EITHER supported source.

    Accepts the raw vendor `<sample>_cell_groups.csv` (cell_id, group, color) or the
    `metadata_<sample>.tsv` written by the `metadata` subcommand (cell id index, one
    label column). Detected from the extension and columns rather than asked for, so
    the caller does not have to care which one is to hand."""
    sep = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    if "group" in df.columns and "cell_id" in df.columns:
        df["cell_id"] = df["cell_id"].astype(str)
        labels = df.set_index("cell_id")["group"]
        kind = "vendor cell_groups"
    else:
        df = df.set_index(df.columns[0])
        df.index = df.index.astype(str)
        col = LABEL_COLUMN if LABEL_COLUMN in df.columns else df.columns[0]
        labels = df[col]
        kind = f"external-annotation TSV (column '{col}')"
    labels = labels.fillna("Unannotated").str.strip()
    print(f"[format_atera] label source: {path}  [{kind}]  "
          f"{len(labels):,} cells, {labels.nunique()} labels")
    return labels


def colour_block(df):
    """`annotation_colors:` YAML using 10x's display colours, alphabetical."""
    if "color" not in df.columns:
        return ""
    pairs = (df.dropna(subset=["color"])
               .groupby("group")["color"]
               .agg(lambda s: s.mode().iat[0])      # one colour per group
               .sort_index())
    lines = ["annotation_colors:"]
    lines += [f"  {g}: '{c}'" for g, c in pairs.items()]
    return "\n".join(lines)


def cmd_metadata(args):
    df = read_cell_groups(args.cell_groups)
    out = pd.DataFrame({LABEL_COLUMN: df["group"]})
    out.index.name = "cell_id"

    os.makedirs(args.outdir, exist_ok=True)
    dest = os.path.join(args.outdir, f"metadata_{args.sample}.tsv")
    out.to_csv(dest, sep="\t")
    print(f"[format_atera] wrote {dest}  ({len(out):,} cells, "
          f"{out[LABEL_COLUMN].nunique()} labels)")

    if args.check_contract:
        try:
            import anndata as ad
            adata = ad.read_h5ad(args.check_contract, backed="r")
            obs_names = set(map(str, adata.obs_names))
            hit = len(obs_names & set(out.index))
            print(f"[format_atera] contract has {len(obs_names):,} cells; "
                  f"{hit:,} ({100 * hit / max(len(obs_names), 1):.1f}%) carry a label")
            if hit == 0:
                print("[format_atera] WARNING: no overlap. The vendor cell ids do not "
                      "match obs_names, so every cell would fall back to "
                      "'Unannotated'. Check how obs_names were set in the contract.")
        except Exception as e:                                     # noqa: BLE001
            print(f"[format_atera] could not check the contract ({e})")

    block = colour_block(df)
    if block:
        print("\n[format_atera] paste into the config:\n")
        print(block)


def cmd_cluster_annotations(args):
    vendor = read_labels(args.labels)

    meta = pd.read_csv(args.pipeline_metadata, sep="\t", index_col=0)
    meta.index = meta.index.astype(str)                 # never let ids become int64
    leiden_col = f"leiden_{str(args.resolution).replace('.', '_')}"
    if leiden_col not in meta.columns:
        avail = [c for c in meta.columns if c.startswith("leiden_")]
        sys.exit(f"[format_atera] {args.pipeline_metadata} has no '{leiden_col}'. "
                 f"Available: {avail or 'none'} — check --resolution.")

    joined = pd.DataFrame({"cluster": meta[leiden_col]}).join(vendor.rename("label"),
                                                              how="left")
    matched = joined["label"].notna().sum()
    print(f"[format_atera] {matched:,}/{len(joined):,} pipeline cells matched a vendor "
          f"label ({100 * matched / max(len(joined), 1):.1f}%)")
    if matched == 0:
        sys.exit("[format_atera] no cell ids matched — nothing to vote on.")

    votes = (joined.dropna(subset=["label"])
                   .groupby("cluster")["label"]
                   .agg(lambda s: s.value_counts().idxmax()))
    purity = (joined.dropna(subset=["label"])
                    .groupby("cluster")["label"]
                    .agg(lambda s: s.value_counts(normalize=True).max()))

    col = f"{args.sample}_{args.resolution}"
    out = pd.DataFrame({col: votes})
    out.index.name = "cluster"
    out = out.sort_index(key=lambda i: pd.to_numeric(i, errors="coerce"))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    out.to_csv(args.out, sep="\t")
    print(f"[format_atera] wrote {args.out}  ({len(out)} clusters)\n")

    print(f"{'cluster':>8}  {'majority label':<40} purity")
    for c in out.index:
        print(f"{c:>8}  {out.loc[c, col]:<40} {purity.get(c, float('nan')):.2f}")
    low = purity[purity < 0.5]
    if len(low):
        print(f"\n[format_atera] NOTE: cluster(s) {list(low.index)} are below 50% "
              "purity — the majority label is a weak summary there. Inspect the "
              "marker plots before trusting them.")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("metadata", help="cell_groups.csv -> external-annotation TSV")
    m.add_argument("cell_groups")
    m.add_argument("--sample", required=True)
    m.add_argument("--outdir", default="demos/atera/external_metadata")
    m.add_argument("--check-contract", metavar="H5AD",
                   help="optional contract h5ad; reports cell-id overlap")
    m.set_defaults(func=cmd_metadata)

    c = sub.add_parser("cluster-annotations",
                       help="label pipeline cells from the vendor file, then "
                            "majority-vote onto the Leiden clusters")
    c.add_argument("--labels", default=None,
                   help="label source: the external-annotation TSV from `metadata` "
                        "(default, resolved from --sample) or the raw cell_groups.csv")
    c.add_argument("--pipeline-metadata", required=True,
                   help="metadata_<sample>.tsv written by preprocess_umap")
    c.add_argument("--sample", required=True)
    c.add_argument("--resolution", default="0.6")
    c.add_argument("--out", default="demos/atera/cluster_annotations.tsv")
    c.set_defaults(func=cmd_cluster_annotations)

    args = p.parse_args()
    if getattr(args, "labels", "sentinel") is None:
        args.labels = os.path.join("demos", "atera", "external_metadata",
                                   f"metadata_{args.sample}.tsv")
        if not os.path.isfile(args.labels):
            sys.exit(f"[format_atera] no label source given and {args.labels} does not "
                     "exist. Run the `metadata` subcommand first, or pass --labels.")
    args.func(args)


if __name__ == "__main__":
    main()
