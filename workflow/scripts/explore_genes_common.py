"""Shared plotting/annotation helpers for the two explore_genes steps
(explore_genes_integrated + explore_genes_sample). Previously duplicated verbatim in
both; kept here as the single definition. Runs in the pseudobulk_aggregate env.
"""
import csv
import logging
import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
from plotting_legends import legend_layout

log = logging.getLogger(__name__)


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
    """Create a temporary legend band composited beneath its dotplots.

    Cell-type and region legends use dedicated rows with adaptive height, so long
    labels cannot overlap. The caller embeds this image in the same output PNG and
    deletes the temporary file; it is never a separate report page.
    """
    try:
        from matplotlib.patches import Patch

        ct_cd = annotation_colors.get(annot_key, {}) if isinstance(annotation_colors, dict) else {}
        reg_cd = region_colors if isinstance(region_colors, dict) else {}
        if not ct_cd and not reg_cd:
            return

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

        groups = []
        if unique_cts:
            groups.append(("Cell type", unique_cts))
        if unique_regs:
            groups.append(("Region", unique_regs))

        layouts = [legend_layout(list(entries), preferred_rows=4)
                   for _, entries in groups]
        heights = [height for _, _, height in layouts]
        fig_leg = plt.figure(figsize=(12, sum(heights)))
        grid = fig_leg.add_gridspec(len(groups), 1, height_ratios=heights, hspace=0.15)
        for row, ((title, entries), (labels, ncol, _)) in enumerate(
            zip(groups, layouts)
        ):
            ax = fig_leg.add_subplot(grid[row, 0])
            ax.set_axis_off()
            patches = [
                Patch(facecolor=color, label=label)
                for label, color in entries.items()
            ]
            ax.legend(
                handles=patches,
                labels=labels,
                title=title,
                loc="upper center",
                fontsize=12,
                title_fontsize=13,
                frameon=True,
                edgecolor="lightgray",
                ncol=ncol,
                handletextpad=0.35,
                columnspacing=0.8,
                labelspacing=0.35,
                borderaxespad=0.0,
            )

        fig_leg.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig_leg)

    except Exception as e:
        log.warning("  Annotation legend image failed: %s", e)
        plt.close("all")


def composite_vertical(image_paths, output_path, dpi=300, pad=30):
    """Stack images vertically with padding, right-aligned, save as PNG."""
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
