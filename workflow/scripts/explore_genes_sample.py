"""
explore_genes_sample.py – Gene/signature exploration (per sample)
=================================================================
Processes a single sample and writes composite PNGs organised by entry:
    gene_exploration/{entry}/Spatial/{sample}_spatial.png
    gene_exploration/{entry}/Dotplots/{sample}_dotplot.png
    gene_exploration/{entry}/Dotplots/{sample}_aucell_dotplot.png  (signatures)

Spatial composite = three panels side by side (expression + cell type + region).
Dotplot composite = three views stacked vertically via PIL (celltype, region,
celltype × region).

Each sample is an independent Snakemake job, parallelised via wildcard.
"""

import gc
import logging
import os
import sys
import tempfile
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image

from composition_barplots import find_niche_column, build_niche_palette
from explore_genes_common import (
    read_tsv_to_dict, classify_entries, apply_annotation_palette, apply_region_palette,
    make_score_adata, annotate_ct_region_dotplot, create_annotation_legend,
    composite_vertical, _compact_legend,
)

# Large cell type × region dotplots can exceed PIL's default pixel limit
Image.MAX_IMAGE_PIXELS = None

sc.settings.set_figure_params(vector_friendly=True)


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
log = logging.getLogger("explore_genes_sample")


# ── Helpers ──────────────────────────────────────────────────────────────────

def generate_spatial_composite(adata, color_col, annot_key, has_regions,
                               library_id, out_path, vmin, vmax, dpi,
                               sample_id, title="expression",
                               niche_col=None, niche_palette=None):
    """Spatial composite: expression/score + cell type + region (+ spatial niche
    when available). The niche panel uses the shared palette and a compacted
    multi-column legend so a high niche count does not distort the composite."""
    has_niche = bool(niche_col and niche_col in adata.obs.columns)
    n_panels = 2 + (1 if has_regions else 0) + (1 if has_niche else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 7),
                              gridspec_kw={"wspace": 0.5})
    fig.suptitle(sample_id, fontsize=16, fontweight="bold", y=1.02)

    # (1) Expression / AUCell
    try:
        sc.pl.spatial(adata, color=color_col, spot_size=20, frameon=False,
                      vmin=vmin, vmax=vmax, cmap="viridis",
                      title=title, library_id=library_id,
                      ax=axes[0], show=False)
    except Exception as e:
        log.warning("    Spatial %s failed: %s", title, e)
        axes[0].set_title(f"{title} (failed)")

    # (2) Cell type
    try:
        sc.pl.spatial(adata, color=annot_key, spot_size=20, frameon=False,
                      title="Cell types", library_id=library_id,
                      legend_fontsize=6, na_in_legend=False,
                      ax=axes[1], show=False)
    except Exception as e:
        log.warning("    Spatial celltype failed: %s", e)
        axes[1].set_title("Cell types (failed)")

    idx = 2
    # (3) Region
    if has_regions:
        try:
            sc.pl.spatial(adata, color="region_annotation", spot_size=20,
                          frameon=False, title="Regions",
                          library_id=library_id,
                          legend_fontsize=6, na_in_legend=False,
                          ax=axes[idx], show=False)
        except Exception as e:
            log.warning("    Spatial region failed: %s", e)
            axes[idx].set_title("Regions (failed)")
        idx += 1

    # (4) Spatial niche — colours consistent with the spatial_niches rule, legend
    #     compacted into columns so a high niche count stays readable.
    if has_niche:
        try:
            if not isinstance(adata.obs[niche_col].dtype, pd.CategoricalDtype):
                adata.obs[niche_col] = adata.obs[niche_col].astype(str).astype("category")
            cats = list(adata.obs[niche_col].cat.categories)
            if niche_palette:
                adata.uns[f"{niche_col}_colors"] = [
                    niche_palette.get(str(c), "#cccccc") for c in cats]
            sc.pl.spatial(adata, color=niche_col, spot_size=20, frameon=False,
                          title="Spatial niches", library_id=library_id,
                          legend_fontsize=5, na_in_legend=False,
                          ax=axes[idx], show=False)
            _compact_legend(axes[idx], title="Niche")
        except Exception as e:
            log.warning("    Spatial niche failed: %s", e)
            axes[idx].set_title("Spatial niches (failed)")
        idx += 1

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def generate_dotplot_composite(adata, var_names, annot_key, has_regions,
                               out_path, prefix, dpi,
                               annotation_colors=None, region_colors=None):
    """Composite PNG of three dotplots stacked vertically via PIL."""
    n_genes = len(var_names)
    fig_w = max(8, n_genes * 1.2 + 4)
    tmp_files = []

    # (1) By cell type
    try:
        n_ct = adata.obs[annot_key].nunique()
        tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
        sc.pl.dotplot(adata, var_names=var_names, groupby=annot_key,
                      standard_scale="var",
                      figsize=(fig_w, max(4, n_ct * 0.4)),
                      title=f"{prefix} – by cell type",
                      show=False)
        plt.savefig(tmp, dpi=dpi, bbox_inches="tight")
        plt.close("all")
        tmp_files.append(tmp)
    except Exception as e:
        log.warning("    Dotplot celltype failed: %s", e)
        plt.close("all")

    # (2) By region
    if has_regions:
        try:
            n_reg = adata.obs["region_annotation"].nunique()
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
            sc.pl.dotplot(adata, var_names=var_names,
                          groupby="region_annotation",
                          standard_scale="var",
                          figsize=(fig_w, max(4, n_reg * 0.5)),
                          title=f"{prefix} – by region",
                          show=False)
            plt.savefig(tmp, dpi=dpi, bbox_inches="tight")
            plt.close("all")
            tmp_files.append(tmp)
        except Exception as e:
            log.warning("    Dotplot region failed: %s", e)
            plt.close("all")

    # (3) By cell type × region
    if has_regions:
        try:
            combined = (adata.obs[annot_key].astype(str) + " | "
                        + adata.obs["region_annotation"].astype(str))
            adata.obs["_ct_region"] = pd.Categorical(combined)
            n_groups = adata.obs["_ct_region"].nunique()
            fig_h = max(6, n_groups * 0.35)
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
            sc.pl.dotplot(adata, var_names=var_names, groupby="_ct_region",
                          standard_scale="var",
                          figsize=(fig_w, fig_h),
                          title=f"{prefix} – by cell type × region",
                          show=False)
            if annotation_colors or region_colors:
                annotate_ct_region_dotplot(annot_key,
                                          annotation_colors or {},
                                          region_colors or {})
            plt.savefig(tmp, dpi=dpi, bbox_inches="tight")
            tmp_files.append(tmp)

            # Standalone annotation legend (before plt.close so gcf() can
            # still read the dotplot's ytick labels)
            if annotation_colors or region_colors:
                tmp_leg = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
                create_annotation_legend(annot_key,
                                         annotation_colors or {},
                                         region_colors or {},
                                         tmp_leg, dpi)
                if os.path.isfile(tmp_leg) and os.path.getsize(tmp_leg) > 0:
                    tmp_files.append(tmp_leg)
            plt.close("all")
        except Exception as e:
            log.warning("    Dotplot celltype×region failed: %s", e)
            plt.close("all")
        finally:
            if "_ct_region" in adata.obs.columns:
                del adata.obs["_ct_region"]

    # Composite
    if tmp_files:
        composite_vertical(tmp_files, out_path, dpi=dpi)
    for f in tmp_files:
        try:
            os.unlink(f)
        except OSError:
            pass


# ── Parameters ───────────────────────────────────────────────────────────────
adata_path        = str(snakemake.input.adata)
ranges_path       = str(snakemake.input.ranges)
queries_path      = str(snakemake.input.queries)
base_dir          = str(snakemake.params.outdir)
sample_id         = str(snakemake.params.sample_id)

ANNOT_KEY         = str(snakemake.params.annot_key)
AUCELL_FRACTION   = float(snakemake.params.aucell_fraction)
DPI               = int(snakemake.params.dpi)
ANNOTATION_COLORS = snakemake.params.annotation_colors
REGION_COLORS     = snakemake.params.region_colors
NICHE_COLUMN      = str(getattr(snakemake.params, "niche_column", "") or "")

try:
    log.info("=" * 70)
    log.info("Gene exploration – sample %s", sample_id)
    log.info("=" * 70)

    # ── 1. Parse queries + ranges ────────────────────────────────────────
    entries_dict = read_tsv_to_dict(queries_path)
    individual_genes, signatures = classify_entries(entries_dict)
    all_entries = {**individual_genes, **signatures}

    ranges_df = pd.read_csv(ranges_path, sep="\t")
    ranges_map = {row["entry"]: (row["vmin"], row["vmax"])
                  for _, row in ranges_df.iterrows()}

    if not all_entries:
        log.warning("No entries — nothing to plot.")
        sys.exit(0)

    # ── 2. Load sample ──────────────────────────────────────────────────
    log.info("Loading %s …", adata_path)
    adata = sc.read_h5ad(adata_path)
    log.info("  %d cells, %d genes", adata.n_obs, adata.n_vars)

    library_id = None
    if "spatial" in adata.uns and len(adata.uns["spatial"]) == 1:
        library_id = list(adata.uns["spatial"].keys())[0]

    apply_annotation_palette(adata, ANNOT_KEY, ANNOTATION_COLORS)
    apply_region_palette(adata, REGION_COLORS)
    if ANNOT_KEY in adata.obs.columns:
        adata.obs[ANNOT_KEY] = adata.obs[ANNOT_KEY].astype("category")

    sample_genes = set(adata.var_names)
    has_regions = ("region_annotation" in adata.obs.columns
                   and adata.obs["region_annotation"].nunique() > 1
                   and not all(adata.obs["region_annotation"] == "Unlabeled"))

    # Spatial niche (if present): shared value-deterministic palette → colours
    # match the spatial_niches rule across the pipeline.
    niche_col = find_niche_column(adata, NICHE_COLUMN or None)
    niche_palette = None
    if niche_col:
        _nv = adata.obs[niche_col].astype(str)
        try:
            _ncats = sorted(_nv.unique(), key=lambda x: int(x))
        except ValueError:
            _ncats = sorted(_nv.unique())
        niche_palette = build_niche_palette(
            _ncats, ANNOTATION_COLORS.get("spatial_niche", {})
            if isinstance(ANNOTATION_COLORS, dict) else {})
        log.info("  spatial niche '%s' (%d niches) → composite panel",
                 niche_col, len(_ncats))

    # ── 3. Validate ──────────────────────────────────────────────────────
    valid_entries = {}
    for name, gene_list in all_entries.items():
        present = [g for g in gene_list if g in sample_genes]
        if present:
            valid_entries[name] = present
    individual_genes = {k: v for k, v in individual_genes.items() if k in valid_entries}
    signatures = {k: v for k, v in signatures.items() if k in valid_entries}

    # ── 4. AUCell ────────────────────────────────────────────────────────
    if signatures:
        log.info("Computing AUCell for %d signatures …", len(signatures))
        import decoupler as dc

        net_rows = []
        for sig_name in signatures:
            for gene in valid_entries[sig_name]:
                net_rows.append({"source": sig_name, "target": gene, "weight": 1.0})
        net_df = pd.DataFrame(net_rows)
        n_up = max(1, int(adata.n_vars * AUCELL_FRACTION))

        try:
            dc.run_aucell(adata, net=net_df, source="source", target="target",
                          n_up=n_up, use_raw=False, verbose=False)
            if "aucell_estimate" in adata.obsm:
                est = adata.obsm["aucell_estimate"]
                for sig_name in signatures:
                    if sig_name in est.columns:
                        adata.obs[f"AUCell_{sig_name}"] = est[sig_name].values
        except Exception as e:
            log.warning("  AUCell failed: %s", e)

    # ── 5. Generate PNGs ─────────────────────────────────────────────────

    for gene_name in individual_genes:
        if gene_name not in ranges_map:
            continue
        log.info("  Gene: %s", gene_name)
        vmin, vmax = ranges_map[gene_name]

        spatial_dir = os.path.join(base_dir, gene_name, "Spatial")
        dotplot_dir = os.path.join(base_dir, gene_name, "Dotplots")
        os.makedirs(spatial_dir, exist_ok=True)
        os.makedirs(dotplot_dir, exist_ok=True)

        # Spatial composite
        generate_spatial_composite(
            adata, gene_name, ANNOT_KEY, has_regions, library_id,
            os.path.join(spatial_dir, f"{sample_id}_spatial.png"),
            vmin, vmax, DPI, sample_id, title="expression",
            niche_col=niche_col, niche_palette=niche_palette)

        # Dotplot composite
        generate_dotplot_composite(
            adata, [gene_name], ANNOT_KEY, has_regions,
            os.path.join(dotplot_dir, f"{sample_id}_dotplot.png"),
            gene_name, DPI, ANNOTATION_COLORS, REGION_COLORS)

    for sig_name in signatures:
        score_col = f"AUCell_{sig_name}"
        if score_col not in adata.obs.columns or sig_name not in ranges_map:
            continue
        log.info("  Signature: %s", sig_name)
        vmin, vmax = ranges_map[sig_name]
        present_genes = valid_entries[sig_name]

        spatial_dir = os.path.join(base_dir, sig_name, "Spatial")
        dotplot_dir = os.path.join(base_dir, sig_name, "Dotplots")
        os.makedirs(spatial_dir, exist_ok=True)
        os.makedirs(dotplot_dir, exist_ok=True)

        # Spatial composite (AUCell)
        generate_spatial_composite(
            adata, score_col, ANNOT_KEY, has_regions, library_id,
            os.path.join(spatial_dir, f"{sample_id}_spatial.png"),
            vmin, vmax, DPI, sample_id, title="AUCell",
            niche_col=niche_col, niche_palette=niche_palette)

        # Dotplot composite – member genes
        generate_dotplot_composite(
            adata, present_genes, ANNOT_KEY, has_regions,
            os.path.join(dotplot_dir, f"{sample_id}_dotplot.png"),
            f"{sig_name} genes", DPI, ANNOTATION_COLORS, REGION_COLORS)

        # Dotplot composite – AUCell score
        score_ad = make_score_adata(adata, score_col)
        generate_dotplot_composite(
            score_ad, [score_col], ANNOT_KEY, has_regions,
            os.path.join(dotplot_dir, f"{sample_id}_aucell_dotplot.png"),
            f"{sig_name} AUCell", DPI, ANNOTATION_COLORS, REGION_COLORS)
        del score_ad

    del adata
    gc.collect()
    log.info("Done – sample %s", sample_id)

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
