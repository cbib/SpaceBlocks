"""
explore_genes_integrated.py – Gene/signature exploration (integrated)
=====================================================================
Phase 1: loads the Harmony-integrated h5ad once, computes AUCell scores
for all signatures, writes expression_ranges.tsv, and produces PNGs
organised as:
    gene_exploration/{entry}/{entry}_genes.tsv
    gene_exploration/{entry}/Integrated/*.png
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
import numpy as np
import pandas as pd
import scanpy as sc

# Rasterise scatter data layers (text/axes stay vector)
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


def make_score_adata(adata, score_col):
    """Minimal AnnData with an obs score as a pseudo-gene for sc.pl.dotplot.
    Drops score_col from obs to avoid 'found in both obs and var' conflict."""
    return sc.AnnData(
        X=np.asarray(adata.obs[score_col].values, dtype=np.float32).reshape(-1, 1),
        obs=adata.obs.drop(columns=[score_col], errors="ignore"),
        var=pd.DataFrame(index=[score_col]),
    )


def annotate_ct_region_dotplot(annot_key, annotation_colors, region_colors):
    """Add cell type and region colour bars + legends to the RIGHT of a
    ct×region dotplot.  Call after sc.pl.dotplot(..., show=False).
    Axes placed beyond the figure boundary are captured by bbox_inches='tight'."""
    try:
        from matplotlib.patches import Patch

        ct_cd = annotation_colors.get(annot_key, {}) if isinstance(annotation_colors, dict) else {}
        reg_cd = region_colors if isinstance(region_colors, dict) else {}
        if not ct_cd and not reg_cd:
            return

        fig = plt.gcf()
        fig.canvas.draw()

        # Find the main axes (ytick labels contain " | ")
        main_ax = None
        for ax in fig.axes:
            labels = [t.get_text() for t in ax.get_yticklabels()]
            if labels and any(" | " in l for l in labels):
                main_ax = ax
                break
        if main_ax is None:
            return

        yticks = main_ax.get_yticks()
        ylabels = [t.get_text() for t in main_ax.get_yticklabels()]
        if not ylabels:
            return

        # Parse labels and collect unique categories
        ct_cols, reg_cols = [], []
        unique_cts, unique_regs = {}, {}
        for label in ylabels:
            parts = label.split(" | ", 1)
            ct = parts[0].strip()
            region = parts[1].strip() if len(parts) > 1 else ""
            ct_c = ct_cd.get(ct, "#cccccc")
            reg_c = reg_cd.get(region, "#cccccc")
            ct_cols.append(ct_c)
            reg_cols.append(reg_c)
            if ct not in unique_cts:
                unique_cts[ct] = ct_c
            if region and region not in unique_regs:
                unique_regs[region] = reg_c

        # Position bars to the RIGHT of all existing content
        max_right = max(ax.get_position().x1 for ax in fig.axes)
        pos = main_ax.get_position()
        ylim = main_ax.get_ylim()
        bar_w = 0.018
        gap = 0.008
        start_x = max_right + gap

        for i, (colors, title) in enumerate([(ct_cols, "Cell type"),
                                              (reg_cols, "Region")]):
            x = start_x + i * (bar_w + gap)
            ann_ax = fig.add_axes([x, pos.y0, bar_w, pos.y1 - pos.y0])
            for yt, col in zip(yticks[:len(colors)], colors):
                ann_ax.barh(yt, 1, height=0.8, color=col, edgecolor="none")
            ann_ax.set_ylim(ylim)
            ann_ax.set_xlim(0, 1)
            ann_ax.set_xticks([])
            ann_ax.set_yticks([])
            ann_ax.set_title(title, fontsize=7)
            for spine in ann_ax.spines.values():
                spine.set_visible(False)

        # Colour legends to the right of the bars
        legend_x = start_x + 2 * (bar_w + gap) + gap
        if unique_cts:
            ct_patches = [Patch(facecolor=c, label=n) for n, c in unique_cts.items()]
            leg1 = fig.legend(handles=ct_patches, title="Cell type",
                              loc="upper left", bbox_to_anchor=(legend_x, pos.y1),
                              fontsize=5, title_fontsize=6,
                              frameon=True, edgecolor="lightgray")
            fig.add_artist(leg1)
        if unique_regs:
            reg_patches = [Patch(facecolor=c, label=n) for n, c in unique_regs.items()]
            fig.legend(handles=reg_patches, title="Region",
                       loc="lower left", bbox_to_anchor=(legend_x, pos.y0),
                       fontsize=5, title_fontsize=6,
                       frameon=True, edgecolor="lightgray")
    except Exception as e:
        log.warning("  Row annotation failed: %s", e)


def generate_dotplots(adata, var_names, annot_key, has_regions, out_dir,
                      prefix, dpi, annotation_colors=None, region_colors=None):
    """Generate three dotplot PNGs: by cell type, by region, cell type × region."""
    n_genes = len(var_names)
    fig_w = max(8, n_genes * 1.2 + 4)

    # (1) By cell type
    try:
        n_ct = adata.obs[annot_key].nunique()
        sc.pl.dotplot(adata, var_names=var_names, groupby=annot_key,
                      standard_scale="var",
                      figsize=(fig_w, max(4, n_ct * 0.4)),
                      title=f"{prefix} – by cell type",
                      show=False)
        plt.savefig(os.path.join(out_dir, f"{prefix}_dotplot_celltype.png"),
                    dpi=dpi, bbox_inches="tight")
        plt.close("all")
    except Exception as e:
        log.warning("  Dotplot celltype failed: %s", e)
        plt.close("all")

    # (2) By region
    if has_regions:
        try:
            n_reg = adata.obs["region_annotation"].nunique()
            sc.pl.dotplot(adata, var_names=var_names,
                          groupby="region_annotation",
                          standard_scale="var",
                          figsize=(fig_w, max(4, n_reg * 0.5)),
                          title=f"{prefix} – by region",
                          show=False)
            plt.savefig(os.path.join(out_dir, f"{prefix}_dotplot_region.png"),
                        dpi=dpi, bbox_inches="tight")
            plt.close("all")
        except Exception as e:
            log.warning("  Dotplot region failed: %s", e)
            plt.close("all")

    # (3) By cell type × region (with row annotations)
    if has_regions:
        try:
            combined = (adata.obs[annot_key].astype(str) + " | "
                        + adata.obs["region_annotation"].astype(str))
            adata.obs["_ct_region"] = pd.Categorical(combined)
            n_groups = adata.obs["_ct_region"].nunique()
            fig_h = max(6, n_groups * 0.35)
            sc.pl.dotplot(adata, var_names=var_names, groupby="_ct_region",
                          standard_scale="var",
                          figsize=(fig_w, fig_h),
                          title=f"{prefix} – by cell type × region",
                          show=False)
            if annotation_colors or region_colors:
                annotate_ct_region_dotplot(annot_key,
                                          annotation_colors or {},
                                          region_colors or {})
            plt.savefig(os.path.join(out_dir,
                        f"{prefix}_dotplot_celltype_region.png"),
                        dpi=dpi, bbox_inches="tight")
            plt.close("all")
        except Exception as e:
            log.warning("  Dotplot celltype×region failed: %s", e)
            plt.close("all")
        finally:
            if "_ct_region" in adata.obs.columns:
                del adata.obs["_ct_region"]


def generate_umap_composite(adata, color_col, annot_key, has_regions, out_path,
                             vmin, vmax, dpi, title="expression"):
    """Single composite PNG with 2–3 UMAPs side by side."""
    n_panels = 2 + (1 if has_regions else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 7))

    try:
        sc.pl.umap(adata, color=color_col, size=2, frameon=False,
                    vmin=vmin, vmax=vmax, cmap="viridis",
                    title=title, ax=axes[0], show=False)
    except Exception as e:
        log.warning("  UMAP %s failed: %s", title, e)

    try:
        sc.pl.umap(adata, color=annot_key, size=2, frameon=False,
                    title="Cell types", ax=axes[1], show=False)
    except Exception as e:
        log.warning("  UMAP celltype failed: %s", e)

    if has_regions:
        try:
            sc.pl.umap(adata, color="region_annotation", size=2, frameon=False,
                        title="Regions", ax=axes[2], show=False)
        except Exception as e:
            log.warning("  UMAP region failed: %s", e)

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ── Parameters ───────────────────────────────────────────────────────────────
integrated_path   = str(snakemake.input.integrated)
queries_path      = str(snakemake.input.queries)
out_ranges        = str(snakemake.output.ranges)
base_dir          = str(snakemake.params.outdir)

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

    # ── 1. Parse queries ─────────────────────────────────────────────────
    entries_dict = read_tsv_to_dict(queries_path)
    individual_genes, signatures = classify_entries(entries_dict)
    all_entries = {**individual_genes, **signatures}
    log.info("Parsed %d entries: %d genes, %d signatures",
             len(all_entries), len(individual_genes), len(signatures))

    if not all_entries:
        log.warning("No entries — writing empty ranges and exiting.")
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
            log.warning("  %s: ALL genes missing — skipping.", name)
    individual_genes = {k: v for k, v in individual_genes.items() if k in valid_entries}
    signatures = {k: v for k, v in signatures.items() if k in valid_entries}
    log.info("  %d valid entries", len(valid_entries))

    # ── 4. AUCell for signatures ─────────────────────────────────────────
    aucell_scores = {}
    if signatures:
        log.info("Computing AUCell for %d signatures …", len(signatures))
        import decoupler as dc

        # Build net DataFrame: source, target, weight
        net_rows = []
        for sig_name in signatures:
            for gene in valid_entries[sig_name]:
                net_rows.append({"source": sig_name, "target": gene, "weight": 1.0})
        net_df = pd.DataFrame(net_rows)
        n_up = max(1, int(adata.n_vars * AUCELL_FRACTION))
        log.info("  n_up = %d (%.1f%% of %d genes)",
                 n_up, AUCELL_FRACTION * 100, adata.n_vars)

        dc.run_aucell(adata, net=net_df, source="source", target="target",
                      n_up=n_up, use_raw=False, verbose=True)

        # Extract scores from obsm
        if "aucell_estimate" in adata.obsm:
            est = adata.obsm["aucell_estimate"]
            for sig_name in signatures:
                if sig_name in est.columns:
                    aucell_scores[sig_name] = est[sig_name].values
                    adata.obs[f"AUCell_{sig_name}"] = est[sig_name].values

    # ── 5. Expression ranges ─────────────────────────────────────────────
    log.info("Computing expression ranges …")
    ranges_rows = []

    # Collect ALL unique genes (standalone + signature members)
    all_unique_genes = set(individual_genes.keys())
    for sig_name in signatures:
        for gene in valid_entries[sig_name]:
            all_unique_genes.add(gene)

    for gene_name in sorted(all_unique_genes):
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
    log.info("  Written → %s", out_ranges)

    # ── 6. Apply palettes ────────────────────────────────────────────────
    apply_annotation_palette(adata, ANNOT_KEY, ANNOTATION_COLORS)
    apply_region_palette(adata, REGION_COLORS)
    if ANNOT_KEY in adata.obs.columns:
        adata.obs[ANNOT_KEY] = adata.obs[ANNOT_KEY].astype("category")

    has_regions = ("region_annotation" in adata.obs.columns
                   and adata.obs["region_annotation"].nunique() > 1
                   and not all(adata.obs["region_annotation"] == "Unlabeled"))
    has_niche = bool(NICHE_COLUMN) and NICHE_COLUMN in adata.obs.columns

    # ── 7. Generate PNGs ─────────────────────────────────────────────────

    for gene_name in individual_genes:
        log.info("Plotting gene: %s", gene_name)
        entry_dir = os.path.join(base_dir, gene_name)
        int_dir = os.path.join(entry_dir, "Integrated")
        os.makedirs(int_dir, exist_ok=True)

        # Provenance TSV
        pd.DataFrame({gene_name: [gene_name]}).to_csv(
            os.path.join(entry_dir, f"{gene_name}_genes.tsv"),
            sep="\t", index=False)

        row = ranges_df[ranges_df["entry"] == gene_name].iloc[0]
        vmin, vmax = row["vmin"], row["vmax"]

        generate_dotplots(adata, [gene_name], ANNOT_KEY, has_regions,
                          int_dir, gene_name, DPI, ANNOTATION_COLORS, REGION_COLORS)
        generate_umap_composite(adata, gene_name, ANNOT_KEY, has_regions,
                                os.path.join(int_dir, f"{gene_name}_UMAPs.png"),
                                vmin, vmax, DPI, "expression")

        if has_niche:
            try:
                n_niche = adata.obs[NICHE_COLUMN].nunique()
                sc.pl.dotplot(adata, var_names=[gene_name],
                              groupby=NICHE_COLUMN, standard_scale="var",
                              figsize=(8, max(4, n_niche * 0.5)),
                              title=f"{gene_name} – by niche", show=False)
                plt.savefig(os.path.join(int_dir,
                            f"{gene_name}_dotplot_niche.png"),
                            dpi=DPI, bbox_inches="tight")
                plt.close("all")
            except Exception as e:
                log.warning("  Niche dotplot failed: %s", e)
                plt.close("all")

    for sig_name in signatures:
        log.info("Plotting signature: %s", sig_name)
        score_col = f"AUCell_{sig_name}"
        if score_col not in adata.obs.columns:
            log.warning("  AUCell score not found — skipping.")
            continue

        entry_dir = os.path.join(base_dir, sig_name)
        int_dir = os.path.join(entry_dir, "Integrated")
        os.makedirs(int_dir, exist_ok=True)

        present_genes = valid_entries[sig_name]
        pd.DataFrame({sig_name: present_genes}).to_csv(
            os.path.join(entry_dir, f"{sig_name}_genes.tsv"),
            sep="\t", index=False)

        row = ranges_df[ranges_df["entry"] == sig_name].iloc[0]
        vmin, vmax = row["vmin"], row["vmax"]

        # Member-gene dotplots
        generate_dotplots(adata, present_genes, ANNOT_KEY, has_regions,
                          int_dir, f"{sig_name}_genes", DPI,
                          ANNOTATION_COLORS, REGION_COLORS)

        # AUCell score dotplots
        score_ad = make_score_adata(adata, score_col)
        generate_dotplots(score_ad, [score_col], ANNOT_KEY, has_regions,
                          int_dir, f"{sig_name}_aucell", DPI,
                          ANNOTATION_COLORS, REGION_COLORS)
        del score_ad

        # UMAPs
        generate_umap_composite(adata, score_col, ANNOT_KEY, has_regions,
                                os.path.join(int_dir, f"{sig_name}_UMAPs.png"),
                                vmin, vmax, DPI, "AUCell")

        if has_niche:
            try:
                score_ad = make_score_adata(adata, score_col)
                n_niche = adata.obs[NICHE_COLUMN].nunique()
                sc.pl.dotplot(score_ad, var_names=[score_col],
                              groupby=NICHE_COLUMN, standard_scale="var",
                              figsize=(8, max(4, n_niche * 0.5)),
                              title=f"{sig_name} AUCell – by niche",
                              show=False)
                plt.savefig(os.path.join(int_dir,
                            f"{sig_name}_aucell_dotplot_niche.png"),
                            dpi=DPI, bbox_inches="tight")
                plt.close("all")
                del score_ad
            except Exception as e:
                log.warning("  Niche dotplot failed: %s", e)
                plt.close("all")

    del adata
    gc.collect()
    log.info("Done.")

except Exception:
    log.error("FAILED:\n%s", traceback.format_exc())
    raise
