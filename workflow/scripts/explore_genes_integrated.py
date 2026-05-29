"""
explore_genes_integrated.py – Gene/signature exploration (integrated)
=====================================================================
Phase 1 of gene exploration.  Loads the Harmony-integrated h5ad once,
computes AUCell scores for all signatures, writes expression_ranges.tsv,
and produces one integrated PDF per entry.

Individual genes  → raw normalised expression for dotplot/violin/UMAP.
Signatures        → AUCell score for violin/UMAP; dotplot of member genes.
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
import matplotlib.gridspec as gridspec
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
log = logging.getLogger("explore_genes_integrated")


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
    """
    Classify entries as 'gene' or 'signature'.

    A column with exactly one gene whose name matches the column header
    is an individual gene.  Everything else is a signature.
    """
    genes, signatures = {}, {}
    for name, gene_list in entries_dict.items():
        if len(gene_list) == 1 and gene_list[0] == name:
            genes[name] = gene_list
        else:
            signatures[name] = gene_list
    return genes, signatures


def apply_annotation_palette(adata, obs_key, annotation_colors):
    """Apply annotation_colors palette to an obs column."""
    if not isinstance(annotation_colors, dict):
        return
    cd = annotation_colors.get(obs_key, {})
    if cd and obs_key in adata.obs.columns:
        adata.obs[obs_key] = adata.obs[obs_key].astype("category")
        cats = adata.obs[obs_key].cat.categories
        adata.uns[f"{obs_key}_colors"] = [cd.get(str(c), "#cccccc") for c in cats]


def apply_region_palette(adata, region_colors):
    """Apply region_colors palette to region_annotation."""
    if region_colors and "region_annotation" in adata.obs.columns:
        adata.obs["region_annotation"] = adata.obs["region_annotation"].astype("category")
        cats = adata.obs["region_annotation"].cat.categories
        adata.uns["region_annotation_colors"] = [
            region_colors.get(str(c), "#cccccc") for c in cats
        ]


# ── Parameters ───────────────────────────────────────────────────────────────
integrated_path   = str(snakemake.input.integrated)
queries_path      = str(snakemake.input.queries)
out_ranges        = str(snakemake.output.ranges)
out_dir           = str(snakemake.output.outdir)

ANNOT_KEY         = str(snakemake.params.annot_key)
AUCELL_FRACTION   = float(snakemake.params.aucell_fraction)
NICHE_COLUMN      = str(snakemake.params.niche_column) if snakemake.params.niche_column else ""
DPI               = int(snakemake.params.dpi)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors

try:
    log.info("=" * 70)
    log.info("Gene exploration – integrated plots")
    log.info("=" * 70)

    os.makedirs(out_dir, exist_ok=True)

    # ── 1. Parse queries ─────────────────────────────────────────────────
    entries_dict = read_tsv_to_dict(queries_path)
    individual_genes, signatures = classify_entries(entries_dict)
    all_entries = {**individual_genes, **signatures}
    log.info("Parsed %d entries: %d individual genes, %d signatures",
             len(all_entries), len(individual_genes), len(signatures))

    if not all_entries:
        log.warning("No entries in query file — writing empty ranges and exiting.")
        pd.DataFrame(columns=["entry", "vmin", "vmax", "type"]).to_csv(
            out_ranges, sep="\t", index=False)
        sys.exit(0)

    # ── 2. Load integrated h5ad ──────────────────────────────────────────
    log.info("Loading integrated h5ad …")
    adata = sc.read_h5ad(integrated_path)
    log.info("  %d cells, %d genes", adata.n_obs, adata.n_vars)

    # ── 3. Validate genes ────────────────────────────────────────────────
    available_genes = set(adata.var_names)
    valid_entries = {}
    for name, gene_list in all_entries.items():
        present = [g for g in gene_list if g in available_genes]
        missing = [g for g in gene_list if g not in available_genes]
        if missing:
            log.warning("  %s: %d/%d genes missing: %s",
                        name, len(missing), len(gene_list), ", ".join(missing))
        if present:
            valid_entries[name] = present
        else:
            log.warning("  %s: ALL genes missing — skipping entirely.", name)
    individual_genes = {k: v for k, v in individual_genes.items() if k in valid_entries}
    signatures = {k: v for k, v in signatures.items() if k in valid_entries}
    log.info("  %d valid entries after gene check", len(valid_entries))

    # ── 4. Compute AUCell for signatures ─────────────────────────────────
    aucell_scores = {}
    if signatures:
        log.info("Computing AUCell scores for %d signatures …", len(signatures))
        import decoupler as dc

        # Build net DataFrame: source, target, weight
        net_rows = []
        for sig_name, gene_list in signatures.items():
            for gene in valid_entries[sig_name]:
                net_rows.append({"source": sig_name, "target": gene, "weight": 1.0})
        net_df = pd.DataFrame(net_rows)

        n_up = max(1, int(adata.n_vars * AUCELL_FRACTION))
        log.info("  AUCell n_up = %d (%.1f%% of %d genes)",
                 n_up, AUCELL_FRACTION * 100, adata.n_vars)

        dc.run_aucell(
            adata, net=net_df, source="source", target="target",
            n_up=n_up, use_raw=False, verbose=True,
        )

        # Extract scores from obsm
        if "aucell_estimate" in adata.obsm:
            est = adata.obsm["aucell_estimate"]
            for sig_name in signatures:
                if sig_name in est.columns:
                    aucell_scores[sig_name] = est[sig_name].values
                    adata.obs[f"AUCell_{sig_name}"] = est[sig_name].values
                    log.info("  %s: AUCell range [%.4f, %.4f]",
                             sig_name,
                             aucell_scores[sig_name].min(),
                             aucell_scores[sig_name].max())

    # ── 5. Expression ranges (p1, p99) ───────────────────────────────────
    log.info("Computing expression ranges …")
    ranges_rows = []
    for gene_name in individual_genes:
        expr = adata[:, gene_name].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        expr = np.asarray(expr).flatten()
        vmin, vmax = float(np.percentile(expr, 1)), float(np.percentile(expr, 99))
        ranges_rows.append({"entry": gene_name, "vmin": vmin, "vmax": vmax, "type": "gene"})

    for sig_name in signatures:
        if sig_name in aucell_scores:
            scores = aucell_scores[sig_name]
            vmin = float(np.percentile(scores, 1))
            vmax = float(np.percentile(scores, 99))
            ranges_rows.append({"entry": sig_name, "vmin": vmin, "vmax": vmax,
                                "type": "signature"})

    ranges_df = pd.DataFrame(ranges_rows)
    ranges_df.to_csv(out_ranges, sep="\t", index=False)
    log.info("  Expression ranges written → %s", out_ranges)

    # ── 6. Apply palettes ────────────────────────────────────────────────
    apply_annotation_palette(adata, ANNOT_KEY, ANNOTATION_COLORS)
    apply_region_palette(adata, REGION_COLORS)

    has_regions = ("region_annotation" in adata.obs.columns
                   and adata.obs["region_annotation"].nunique() > 1
                   and not all(adata.obs["region_annotation"] == "Unlabeled"))

    has_niche = bool(NICHE_COLUMN) and NICHE_COLUMN in adata.obs.columns

    # Ensure annot_key is categorical
    if ANNOT_KEY in adata.obs.columns:
        adata.obs[ANNOT_KEY] = adata.obs[ANNOT_KEY].astype("category")

    # ── 7. Generate PDFs ─────────────────────────────────────────────────

    # ---- Individual genes ----
    for gene_name in individual_genes:
        log.info("Plotting gene: %s", gene_name)
        gene_range = ranges_df[ranges_df["entry"] == gene_name].iloc[0]
        vmin, vmax = gene_range["vmin"], gene_range["vmax"]

        try:
            pdf_path = os.path.join(out_dir, f"{gene_name}_integrated.pdf")
            with PdfPages(pdf_path) as pdf:

                # ── Page 1: Dotplot + Violin by cell type ────────────
                fig = plt.figure(figsize=(18, 8))
                fig.suptitle(f"{gene_name} – Cell type overview", fontsize=16,
                             fontweight="bold", y=0.98)
                gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.4)

                ax1 = fig.add_subplot(gs[0, 0])
                try:
                    sc.pl.dotplot(adata, var_names=[gene_name], groupby=ANNOT_KEY,
                                  standard_scale="var", swap_axes=True,
                                  ax=ax1, show=False)
                except Exception as e:
                    log.warning("  Dotplot failed: %s", e)
                    ax1.text(0.5, 0.5, f"Dotplot failed:\n{e}",
                             ha="center", va="center", transform=ax1.transAxes)

                ax2 = fig.add_subplot(gs[0, 1])
                try:
                    sc.pl.violin(adata, keys=gene_name, groupby=ANNOT_KEY,
                                 rotation=90, ax=ax2, show=False)
                except Exception as e:
                    log.warning("  Violin failed: %s", e)
                    ax2.text(0.5, 0.5, f"Violin failed:\n{e}",
                             ha="center", va="center", transform=ax2.transAxes)

                pdf.savefig(fig, dpi=DPI, bbox_inches="tight")
                plt.close(fig)

                # ── Page 2: UMAP expression + cell type + region ─────
                n_cols = 2 + (1 if has_regions else 0)
                fig, axes = plt.subplots(1, n_cols, figsize=(8 * n_cols, 7))
                if n_cols == 1:
                    axes = [axes]
                fig.suptitle(f"{gene_name} – UMAP", fontsize=16,
                             fontweight="bold", y=0.98)

                try:
                    sc.pl.umap(adata, color=gene_name, size=2, frameon=False,
                               vmin=vmin, vmax=vmax, cmap="viridis",
                               title=f"{gene_name} expression",
                               ax=axes[0], show=False)
                except Exception as e:
                    log.warning("  UMAP expression failed: %s", e)

                try:
                    sc.pl.umap(adata, color=ANNOT_KEY, size=2, frameon=False,
                               title="Cell types", ax=axes[1], show=False)
                except Exception as e:
                    log.warning("  UMAP cell type failed: %s", e)

                if has_regions:
                    try:
                        sc.pl.umap(adata, color="region_annotation", size=2,
                                   frameon=False, title="Regions",
                                   ax=axes[2], show=False)
                    except Exception as e:
                        log.warning("  UMAP region failed: %s", e)

                pdf.savefig(fig, dpi=DPI, bbox_inches="tight")
                plt.close(fig)

                # ── Page 3: Regional view ────────────────────────────
                if has_regions:
                    fig = plt.figure(figsize=(18, 8))
                    fig.suptitle(f"{gene_name} – Regional view", fontsize=16,
                                 fontweight="bold", y=0.98)
                    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.4)

                    ax_v = fig.add_subplot(gs[0, 0])
                    try:
                        sc.pl.violin(adata, keys=gene_name,
                                     groupby="region_annotation",
                                     rotation=45, ax=ax_v, show=False)
                        ax_v.set_title(f"{gene_name} by region")
                    except Exception as e:
                        log.warning("  Region violin failed: %s", e)

                    ax_d = fig.add_subplot(gs[0, 1])
                    try:
                        sc.pl.dotplot(adata, var_names=[gene_name],
                                      groupby="region_annotation",
                                      standard_scale="var", swap_axes=True,
                                      ax=ax_d, show=False)
                    except Exception as e:
                        log.warning("  Region dotplot failed: %s", e)

                    pdf.savefig(fig, dpi=DPI, bbox_inches="tight")
                    plt.close(fig)

                # ── Page 4: Niche view (conditional) ─────────────────
                if has_niche:
                    fig, ax = plt.subplots(figsize=(12, 6))
                    fig.suptitle(f"{gene_name} – Niche view", fontsize=16,
                                 fontweight="bold")
                    try:
                        sc.pl.violin(adata, keys=gene_name,
                                     groupby=NICHE_COLUMN,
                                     rotation=45, ax=ax, show=False)
                    except Exception as e:
                        log.warning("  Niche violin failed: %s", e)

                    pdf.savefig(fig, dpi=DPI, bbox_inches="tight")
                    plt.close(fig)

            log.info("  Saved → %s", pdf_path)

        except Exception as e:
            log.warning("FAILED gene %s: %s\n%s", gene_name, e,
                        traceback.format_exc())

    # ---- Signatures ----
    for sig_name, sig_genes in signatures.items():
        log.info("Plotting signature: %s (%d genes)", sig_name, len(sig_genes))
        score_col = f"AUCell_{sig_name}"
        if score_col not in adata.obs.columns:
            log.warning("  AUCell score not found — skipping.")
            continue

        sig_range = ranges_df[ranges_df["entry"] == sig_name].iloc[0]
        vmin, vmax = sig_range["vmin"], sig_range["vmax"]
        present_genes = valid_entries[sig_name]

        try:
            pdf_path = os.path.join(out_dir, f"{sig_name}_integrated.pdf")
            with PdfPages(pdf_path) as pdf:

                # ── Page 1: Dotplot of member genes + AUCell violin ──
                fig = plt.figure(figsize=(18, 8))
                fig.suptitle(f"{sig_name} – Cell type overview", fontsize=16,
                             fontweight="bold", y=0.98)
                gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.4)

                ax1 = fig.add_subplot(gs[0, 0])
                try:
                    sc.pl.dotplot(adata, var_names=present_genes,
                                  groupby=ANNOT_KEY,
                                  standard_scale="var", swap_axes=True,
                                  title=f"{sig_name} genes",
                                  ax=ax1, show=False)
                except Exception as e:
                    log.warning("  Dotplot failed: %s", e)
                    ax1.text(0.5, 0.5, f"Dotplot failed:\n{e}",
                             ha="center", va="center", transform=ax1.transAxes)

                ax2 = fig.add_subplot(gs[0, 1])
                try:
                    sc.pl.violin(adata, keys=score_col, groupby=ANNOT_KEY,
                                 rotation=90, ax=ax2, show=False)
                    ax2.set_title(f"AUCell score by cell type")
                    ax2.set_ylabel("AUCell score")
                except Exception as e:
                    log.warning("  AUCell violin failed: %s", e)
                    ax2.text(0.5, 0.5, f"Violin failed:\n{e}",
                             ha="center", va="center", transform=ax2.transAxes)

                pdf.savefig(fig, dpi=DPI, bbox_inches="tight")
                plt.close(fig)

                # ── Page 2: UMAP AUCell score + cell type + region ───
                n_cols = 2 + (1 if has_regions else 0)
                fig, axes = plt.subplots(1, n_cols, figsize=(8 * n_cols, 7))
                if n_cols == 1:
                    axes = [axes]
                fig.suptitle(f"{sig_name} – UMAP", fontsize=16,
                             fontweight="bold", y=0.98)

                try:
                    sc.pl.umap(adata, color=score_col, size=2, frameon=False,
                               vmin=vmin, vmax=vmax, cmap="viridis",
                               title=f"{sig_name} AUCell",
                               ax=axes[0], show=False)
                except Exception as e:
                    log.warning("  UMAP AUCell failed: %s", e)

                try:
                    sc.pl.umap(adata, color=ANNOT_KEY, size=2, frameon=False,
                               title="Cell types", ax=axes[1], show=False)
                except Exception as e:
                    log.warning("  UMAP cell type failed: %s", e)

                if has_regions:
                    try:
                        sc.pl.umap(adata, color="region_annotation", size=2,
                                   frameon=False, title="Regions",
                                   ax=axes[2], show=False)
                    except Exception as e:
                        log.warning("  UMAP region failed: %s", e)

                pdf.savefig(fig, dpi=DPI, bbox_inches="tight")
                plt.close(fig)

                # ── Page 3: Regional view ────────────────────────────
                if has_regions:
                    fig = plt.figure(figsize=(18, 8))
                    fig.suptitle(f"{sig_name} – Regional view", fontsize=16,
                                 fontweight="bold", y=0.98)
                    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.4)

                    ax_v = fig.add_subplot(gs[0, 0])
                    try:
                        sc.pl.violin(adata, keys=score_col,
                                     groupby="region_annotation",
                                     rotation=45, ax=ax_v, show=False)
                        ax_v.set_title(f"AUCell score by region")
                        ax_v.set_ylabel("AUCell score")
                    except Exception as e:
                        log.warning("  Region violin failed: %s", e)

                    ax_d = fig.add_subplot(gs[0, 1])
                    try:
                        sc.pl.dotplot(adata, var_names=present_genes,
                                      groupby="region_annotation",
                                      standard_scale="var", swap_axes=True,
                                      title=f"{sig_name} genes by region",
                                      ax=ax_d, show=False)
                    except Exception as e:
                        log.warning("  Region dotplot failed: %s", e)

                    pdf.savefig(fig, dpi=DPI, bbox_inches="tight")
                    plt.close(fig)

                # ── Page 4: Niche view (conditional) ─────────────────
                if has_niche:
                    fig, ax = plt.subplots(figsize=(12, 6))
                    fig.suptitle(f"{sig_name} – Niche view", fontsize=16,
                                 fontweight="bold")
                    try:
                        sc.pl.violin(adata, keys=score_col,
                                     groupby=NICHE_COLUMN,
                                     rotation=45, ax=ax, show=False)
                        ax.set_ylabel("AUCell score")
                    except Exception as e:
                        log.warning("  Niche violin failed: %s", e)

                    pdf.savefig(fig, dpi=DPI, bbox_inches="tight")
                    plt.close(fig)

            log.info("  Saved → %s", pdf_path)

        except Exception as e:
            log.warning("FAILED signature %s: %s\n%s", sig_name, e,
                        traceback.format_exc())

    del adata
    gc.collect()
    log.info("Done.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
