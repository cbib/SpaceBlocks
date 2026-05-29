"""
explore_genes_spatial.py – Gene/signature exploration (spatial)
===============================================================
Phase 2 of gene exploration.  Loads each annotated sample sequentially
(one at a time), computes per-sample AUCell scores for signatures, and
appends one spatial page per sample into the corresponding entry's PDF.

Memory is bounded to a single sample at any time.  All PdfPages handles
are lightweight (they write to disk immediately on savefig).
"""

import csv
import gc
import logging
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scanpy as sc


# ── Logging ──────────────────────────────────────────────────────────────────
log_handlers = [logging.StreamHandler(sys.stderr)]
if hasattr(snakemake, "log"):
    if snakemake.log.out:
        Path(snakemake.log.out).parent.mkdir(parents=True, exist_ok=True)
        log_handlers.append(logging.FileHandler(snakemake.log.out, mode="w"))
    if snakemake.log.err:
        Path(snakemake.log.err).parent.mkdir(parents=True, exist_ok=True)
        sys.stderr = open(snakemake.log.err, "w")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=log_handlers,
)
log = logging.getLogger("explore_genes_spatial")


# ── Helpers ──────────────────────────────────────────────────────────────────

def read_tsv_to_dict(tsv_path):
    """Read TSV where columns = entries and rows = genes (pipeline format)."""
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = {field: [] for field in reader.fieldnames}
        for row in reader:
            for key, val in row.items():
                if val.strip():
                    columns[key].append(val.strip())
    return columns


def classify_entries(entries_dict):
    """Single-gene column whose name matches the gene → gene; else → signature."""
    genes, signatures = {}, {}
    for name, gene_list in entries_dict.items():
        if len(gene_list) == 1 and gene_list[0] == name:
            genes[name] = gene_list
        else:
            signatures[name] = gene_list
    return genes, signatures


def apply_annotation_palette(adata, obs_key, annotation_colors):
    if not isinstance(annotation_colors, dict):
        return
    cd = annotation_colors.get(obs_key, {})
    if cd and obs_key in adata.obs.columns:
        adata.obs[obs_key] = adata.obs[obs_key].astype("category")
        cats = adata.obs[obs_key].cat.categories
        adata.uns[f"{obs_key}_colors"] = [cd.get(str(c), "#cccccc") for c in cats]


def apply_region_palette(adata, region_colors):
    if region_colors and "region_annotation" in adata.obs.columns:
        adata.obs["region_annotation"] = adata.obs["region_annotation"].astype("category")
        cats = adata.obs["region_annotation"].cat.categories
        adata.uns["region_annotation_colors"] = [
            region_colors.get(str(c), "#cccccc") for c in cats
        ]


# ── Parameters ───────────────────────────────────────────────────────────────
annotated_paths   = [str(p) for p in snakemake.input.annotated]
ranges_path       = str(snakemake.input.ranges)
queries_path      = str(snakemake.input.queries)
sample_ids        = list(snakemake.params.sample_ids)

ANNOT_KEY         = str(snakemake.params.annot_key)
AUCELL_FRACTION   = float(snakemake.params.aucell_fraction)
DPI               = int(snakemake.params.dpi)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors

# Output directory for spatial PDFs (sibling of .done sentinel)
out_dir = os.path.dirname(str(snakemake.output.done))

try:
    log.info("=" * 70)
    log.info("Gene exploration – spatial plots for %d samples", len(sample_ids))
    log.info("=" * 70)

    os.makedirs(out_dir, exist_ok=True)

    # ── 1. Parse queries + expression ranges ─────────────────────────────
    entries_dict = read_tsv_to_dict(queries_path)
    individual_genes, signatures = classify_entries(entries_dict)
    all_entries = {**individual_genes, **signatures}

    ranges_df = pd.read_csv(ranges_path, sep="\t")
    ranges_map = {row["entry"]: (row["vmin"], row["vmax"])
                  for _, row in ranges_df.iterrows()}

    log.info("  %d entries, %d expression ranges loaded",
             len(all_entries), len(ranges_map))

    if not all_entries:
        log.warning("No entries — nothing to plot.")
        sys.exit(0)

    # ── 2. Open PdfPages handles ─────────────────────────────────────────
    pdf_handles = {}
    for name in all_entries:
        if name in ranges_map:  # only entries with valid ranges
            pdf_path = os.path.join(out_dir, f"{name}_spatial.pdf")
            pdf_handles[name] = PdfPages(pdf_path)

    log.info("  Opened %d PDF handles", len(pdf_handles))

    # Prepare AUCell net DataFrame if there are signatures
    net_df = None
    if signatures:
        import decoupler as dc
        net_rows = []
        for sig_name, gene_list in signatures.items():
            if sig_name in pdf_handles:
                for gene in gene_list:
                    net_rows.append({"source": sig_name, "target": gene, "weight": 1.0})
        if net_rows:
            net_df = pd.DataFrame(net_rows)

    # ── 3. Process samples sequentially ──────────────────────────────────
    for adata_path, sample_id in zip(annotated_paths, sample_ids):
        log.info("Processing sample: %s", sample_id)
        try:
            adata = sc.read_h5ad(adata_path)
            log.info("  Loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)

            # Library ID for spatial plots
            library_id = None
            if "spatial" in adata.uns and len(adata.uns["spatial"]) == 1:
                library_id = list(adata.uns["spatial"].keys())[0]

            # Apply palettes
            apply_annotation_palette(adata, ANNOT_KEY, ANNOTATION_COLORS)
            apply_region_palette(adata, REGION_COLORS)

            # Ensure annot_key is categorical
            if ANNOT_KEY in adata.obs.columns:
                adata.obs[ANNOT_KEY] = adata.obs[ANNOT_KEY].astype("category")

            # Available genes in this sample
            sample_genes = set(adata.var_names)

            # Compute AUCell scores for this sample
            if net_df is not None:
                # Filter net to genes present in this sample
                net_sample = net_df[net_df["target"].isin(sample_genes)]
                if len(net_sample) > 0:
                    n_up = max(1, int(adata.n_vars * AUCELL_FRACTION))
                    try:
                        dc.run_aucell(
                            adata, net=net_sample, source="source", target="target",
                            n_up=n_up, use_raw=False, verbose=False,
                        )
                        if "aucell_estimate" in adata.obsm:
                            est = adata.obsm["aucell_estimate"]
                            for sig_name in signatures:
                                if sig_name in est.columns:
                                    adata.obs[f"AUCell_{sig_name}"] = est[sig_name].values
                    except Exception as e:
                        log.warning("  AUCell failed for %s: %s", sample_id, e)

            # ── Plot each entry ──────────────────────────────────────
            for name in list(pdf_handles.keys()):
                try:
                    vmin, vmax = ranges_map.get(name, (None, None))
                    is_signature = name in signatures

                    if is_signature:
                        color_col = f"AUCell_{name}"
                        if color_col not in adata.obs.columns:
                            continue
                        title_left = f"{name} AUCell"
                    else:
                        color_col = name
                        if name not in sample_genes:
                            continue
                        title_left = f"{name} expression"

                    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
                    fig.suptitle(f"{sample_id}", fontsize=16,
                                 fontweight="bold", y=0.98)

                    # Left: spatial expression / AUCell score
                    try:
                        sc.pl.spatial(
                            adata, color=color_col, spot_size=20, frameon=False,
                            vmin=vmin, vmax=vmax, cmap="viridis",
                            title=title_left,
                            library_id=library_id, ax=axes[0], show=False,
                        )
                    except Exception as e:
                        log.warning("    Spatial %s failed: %s", name, e)
                        axes[0].text(0.5, 0.5, f"Failed: {e}",
                                     ha="center", va="center",
                                     transform=axes[0].transAxes, fontsize=8)
                        axes[0].set_title(title_left)

                    # Right: spatial cell type annotation
                    try:
                        sc.pl.spatial(
                            adata, color=ANNOT_KEY, spot_size=20, frameon=False,
                            title="Cell types",
                            library_id=library_id, ax=axes[1], show=False,
                        )
                    except Exception as e:
                        log.warning("    Spatial cell type failed: %s", e)
                        axes[1].text(0.5, 0.5, f"Failed: {e}",
                                     ha="center", va="center",
                                     transform=axes[1].transAxes, fontsize=8)
                        axes[1].set_title("Cell types")

                    pdf_handles[name].savefig(fig, dpi=DPI, bbox_inches="tight")
                    plt.close(fig)

                except Exception as e:
                    log.warning("  Entry %s failed for %s: %s", name, sample_id, e)

        except Exception as e:
            log.warning("FAILED sample %s: %s\n%s", sample_id, e,
                        traceback.format_exc())
        finally:
            if "adata" in dir():
                del adata
            gc.collect()

    # ── 4. Close all PDFs ────────────────────────────────────────────────
    for name, handle in pdf_handles.items():
        handle.close()
        log.info("  Saved → %s", os.path.join(out_dir, f"{name}_spatial.pdf"))

    log.info("Done.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    # Close any open handles before re-raising
    for handle in pdf_handles.values():
        try:
            handle.close()
        except Exception:
            pass
    raise
