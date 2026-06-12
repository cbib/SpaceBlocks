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

import csv
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


def build_entry_score_adata(adata, individual_genes, signatures):
    """AnnData with one var per ENTRY: log-norm expression for individual genes,
    the AUCell score for signatures (signature member genes are NOT shown
    individually). Returns None if nothing usable."""
    names, cols = [], []
    for g in individual_genes:
        if g in adata.var_names:
            x = adata[:, g].X
            x = x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
            names.append(g)
            cols.append(np.asarray(x, dtype=np.float32))
    for s in signatures:
        score_col = f"AUCell_{s}"
        if score_col in adata.obs.columns:
            names.append(s)
            cols.append(np.asarray(adata.obs[score_col].values, dtype=np.float32))
    if not names:
        return None
    X = np.column_stack(cols).astype(np.float32)
    obs = adata.obs.drop(columns=[n for n in names if n in adata.obs.columns],
                         errors="ignore").copy()
    return sc.AnnData(X=X, obs=obs, var=pd.DataFrame(index=names))


def generate_score_dotplot_composite(score_ad, annot_key, has_regions,
                                     sample_id, out_path, dpi):
    """Composite of three dotplots (entries on X) for one sample:
    by patient, by patient × area, by patient × cell type. Stacked via PIL."""
    score_ad = score_ad.copy()
    score_ad.obs["_patient"] = str(sample_id)
    groupbys = [("_patient", "by patient")]
    if has_regions and "region_annotation" in score_ad.obs.columns:
        score_ad.obs["_patient_region"] = (
            str(sample_id) + " | " + score_ad.obs["region_annotation"].astype(str))
        groupbys.append(("_patient_region", "by patient \u00d7 area"))
    if annot_key in score_ad.obs.columns:
        score_ad.obs["_patient_ct"] = (
            str(sample_id) + " | " + score_ad.obs[annot_key].astype(str))
        groupbys.append(("_patient_ct", "by patient \u00d7 cell type"))

    var_names = list(score_ad.var_names)
    fig_w = max(8, len(var_names) * 0.45 + 4)
    tmp_files = []
    for gb, title in groupbys:
        try:
            n_groups = score_ad.obs[gb].nunique()
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
            sc.pl.dotplot(score_ad, var_names=var_names, groupby=gb,
                          standard_scale="var",
                          figsize=(fig_w, max(3, n_groups * 0.4)),
                          title=f"{sample_id} \u2013 {title}", show=False)
            plt.savefig(tmp, dpi=dpi, bbox_inches="tight")
            plt.close("all")
            tmp_files.append(tmp)
        except Exception as e:
            log.warning("    Score dotplot (%s) failed: %s", title, e)
            plt.close("all")
    if tmp_files:
        composite_vertical(tmp_files, out_path, dpi=dpi)
    for f in tmp_files:
        try:
            os.unlink(f)
        except OSError:
            pass


def composite_vertical(image_paths, output_path, dpi=300, pad=30):
    """Stack images vertically with padding, save as PNG."""
    imgs = []
    for p in image_paths:
        if os.path.isfile(p):
            imgs.append(Image.open(p))
    if not imgs:
        return

    max_w = max(img.width for img in imgs)
    total_h = sum(img.height for img in imgs) + pad * (len(imgs) - 1)
    composite = Image.new("RGB", (max_w, total_h), (255, 255, 255))

    y = 0
    for img in imgs:
        if img.mode == "RGBA":
            bg = Image.new("RGB", img.size, (255, 255, 255))
            bg.paste(img, mask=img.split()[3])
            img = bg
        composite.paste(img, (max_w - img.width, y))  # right-align
        y += img.height + pad

    composite.save(output_path, dpi=(dpi, dpi))
    for img in imgs:
        img.close()


def generate_spatial_composite(adata, color_col, annot_key, has_regions,
                               library_id, out_path, vmin, vmax, dpi,
                               sample_id, title="expression"):
    """Three-panel spatial composite: expression/score + cell type + region."""
    n_panels = 2 + (1 if has_regions else 0)
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

    # (3) Region
    if has_regions:
        try:
            sc.pl.spatial(adata, color="region_annotation", spot_size=20,
                          frameon=False, title="Regions",
                          library_id=library_id,
                          legend_fontsize=6, na_in_legend=False,
                          ax=axes[2], show=False)
        except Exception as e:
            log.warning("    Spatial region failed: %s", e)
            axes[2].set_title("Regions (failed)")

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def annotate_ct_region_dotplot(annot_key, annotation_colors, region_colors):
    """Add cell type and region colour bars BETWEEN the row labels and the
    dot grid, with colour legends below.  tick_params(pad) creates space
    between labels and axes edge; bars are placed in that gap."""
    try:
        from matplotlib.patches import Patch

        ct_cd = annotation_colors.get(annot_key, {}) if isinstance(annotation_colors, dict) else {}
        reg_cd = region_colors if isinstance(region_colors, dict) else {}
        if not ct_cd and not reg_cd:
            return

        fig = plt.gcf()

        # Find main axes (ytick labels contain " | ")
        main_ax = None
        for ax in fig.axes:
            fig.canvas.draw()
            labels = [t.get_text() for t in ax.get_yticklabels()]
            if labels and any(" | " in l for l in labels):
                main_ax = ax
                break
        if main_ax is None:
            return

        # Widen the gap between labels and axes edge
        main_ax.tick_params(axis="y", pad=30)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

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

        pos = main_ax.get_position()
        ylim = main_ax.get_ylim()

        # Find where labels end (rightmost edge in figure coords)
        max_label_right = 0
        for label_artist in main_ax.get_yticklabels():
            bbox = label_artist.get_window_extent(renderer)
            fig_bbox = bbox.transformed(fig.transFigure.inverted())
            max_label_right = max(max_label_right, fig_bbox.x1)

        # Fit bars in the gap: [labels] gap [CT][Region] gap [dots]
        available = pos.x0 - max_label_right
        bar_w = min(0.020, available * 0.38)
        gap = min(0.006, available * 0.08)

        ct_x = max_label_right + gap
        reg_x = ct_x + bar_w + gap

        for x, colors, title in [(ct_x, ct_cols, "Cell type"),
                                  (reg_x, reg_cols, "Region")]:
            ann_ax = fig.add_axes([x, pos.y0, bar_w, pos.y1 - pos.y0])
            for yt, col in zip(yticks[:len(colors)], colors):
                ann_ax.barh(yt, 1, height=0.8, color=col, edgecolor="none")
            ann_ax.set_ylim(ylim)
            ann_ax.set_xlim(0, 1)
            ann_ax.set_xticks([])
            ann_ax.set_yticks([])
            ann_ax.set_title(title, fontsize=8, rotation=90, ha="left", va="bottom")
            for spine in ann_ax.spines.values():
                spine.set_visible(False)

    except Exception as e:
        log.warning("  Row annotation failed: %s", e)


def create_annotation_legend(annot_key, annotation_colors, region_colors,
                             out_path, dpi):
    """Create a standalone legend image for cell type + region colour bars.
    Legends are stacked vertically so they never overlap."""
    try:
        from matplotlib.patches import Patch

        ct_cd = annotation_colors.get(annot_key, {}) if isinstance(annotation_colors, dict) else {}
        reg_cd = region_colors if isinstance(region_colors, dict) else {}
        if not ct_cd and not reg_cd:
            return

        # Collect unique categories that actually appear in the last dotplot
        fig_prev = plt.gcf()
        main_ax = None
        for ax in fig_prev.axes:
            labels = [t.get_text() for t in ax.get_yticklabels()]
            if labels and any(" | " in l for l in labels):
                main_ax = ax
                break

        unique_cts, unique_regs = {}, {}
        if main_ax is not None:
            for label in [t.get_text() for t in main_ax.get_yticklabels()]:
                parts = label.split(" | ", 1)
                ct = parts[0].strip()
                region = parts[1].strip() if len(parts) > 1 else ""
                if ct and ct not in unique_cts:
                    unique_cts[ct] = ct_cd.get(ct, "#cccccc")
                if region and region not in unique_regs:
                    unique_regs[region] = reg_cd.get(region, "#cccccc")

        if not unique_cts and not unique_regs:
            return

        # Two axes stacked vertically — one legend per row, no overlap
        n_rows = (1 if unique_cts else 0) + (1 if unique_regs else 0)
        fig_leg, axes_leg = plt.subplots(n_rows, 1, figsize=(12, 1.5 * n_rows),
                                          squeeze=False)
        row = 0
        if unique_cts:
            ax = axes_leg[row, 0]
            ax.set_axis_off()
            ct_patches = [Patch(facecolor=c, label=n) for n, c in unique_cts.items()]
            ax.legend(handles=ct_patches, title="Cell type",
                      loc="center", fontsize=12, title_fontsize=13,
                      frameon=True, edgecolor="lightgray",
                      ncol=max(1, len(unique_cts) // 4 + 1))
            row += 1

        if unique_regs:
            ax = axes_leg[row, 0]
            ax.set_axis_off()
            reg_patches = [Patch(facecolor=c, label=n) for n, c in unique_regs.items()]
            ax.legend(handles=reg_patches, title="Region",
                      loc="center", fontsize=12, title_fontsize=13,
                      frameon=True, edgecolor="lightgray",
                      ncol=max(1, len(unique_regs) // 3 + 1))

        fig_leg.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig_leg)

    except Exception as e:
        log.warning("  Annotation legend image failed: %s", e)
        plt.close("all")


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

    # Per-sample score/expression composite (all entries together): one PNG
    # with three dotplots — by patient, by patient × area, by patient × cell type.
    try:
        score_ad_all = build_entry_score_adata(adata, individual_genes, signatures)
        if score_ad_all is not None:
            persample_dir = os.path.join(base_dir, "PerSampleScores")
            os.makedirs(persample_dir, exist_ok=True)
            generate_score_dotplot_composite(
                score_ad_all, ANNOT_KEY, has_regions, sample_id,
                os.path.join(persample_dir, f"{sample_id}_score_dotplots.png"), DPI)
            del score_ad_all
    except Exception as e:
        log.warning("  Per-sample score dotplot composite failed: %s", e)

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
            vmin, vmax, DPI, sample_id, title="expression")

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
            vmin, vmax, DPI, sample_id, title="AUCell")

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
