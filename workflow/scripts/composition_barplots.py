"""
composition_barplots.py — shared stacked composition barplots
=============================================================
Imported by annotate_cells, integrate_samples, subcluster and sample_report.

Every plot produced here guarantees the four properties requested:

  1. absolute (counts) AND relative (%) variants
  2. a thin black outline around every stack segment
  3. colours pulled from the config palette (annotation_colors / region_colors)
  4. legend order top->bottom matches the visual stack order top->bottom.

On (4): pandas/matplotlib stack the FIRST column at the BOTTOM and list it at
the TOP of the legend, so by default the legend is the *reverse* of the visible
stack. This module fixes that by (a) plotting the columns reversed so the first
category in `order` becomes the TOP segment, and (b) building the legend handles
in `order` (first = top). Stack-top and legend-top therefore coincide.

Place this file in the same directory as the other rule scripts
(e.g. workflow/scripts/). Import it from a Snakemake script with:

    import os, sys
    try:
        _here = os.path.dirname(os.path.abspath(__file__))
    except NameError:                      # very old Snakemake
        _here = os.getcwd()
    sys.path.insert(0, _here)
    from composition_barplots import (
        save_stacked_composition, draw_stacked_composition,
        composition_pair, find_niche_column,
    )
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

DEFAULT_COLOR = "#cccccc"
EDGE_COLOR = "black"
EDGE_LW = 0.4


# ── colour resolution ───────────────────────────────────────────────────────

def _default_palette(n):
    """Return a list of n visually-distinct hex colours.

    Cycles tab20 → tab20b → tab20c (60 distinct), then falls back to evenly
    spaced HSV so an arbitrarily long category list never collapses to grey.
    """
    base = []
    for name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(name)
        base += [matplotlib.colors.to_hex(cmap(i)) for i in range(cmap.N)]
    if n > len(base):
        hsv = plt.get_cmap("hsv")
        n_extra = n - len(base)
        base += [matplotlib.colors.to_hex(hsv(i / max(n_extra, 1)))
                 for i in range(n_extra)]
    return base[:n]


def _resolve_palette(categories, color_map):
    """Map every category to a colour: config palette wins, the remainder get
    distinct default colours (so an empty/partial color_map is never all-grey).
    """
    color_map = color_map if isinstance(color_map, dict) else {}
    missing = [c for c in categories if str(c) not in color_map]
    defaults = dict(zip((str(c) for c in missing), _default_palette(len(missing))))
    return {str(c): color_map.get(str(c), defaults.get(str(c), DEFAULT_COLOR))
            for c in categories}


# ── ordering ────────────────────────────────────────────────────────────────

def _ordered_columns(ct, color_map):
    """Decide the top->bottom category order.

    If a palette dict is supplied its key order wins (this is what makes the
    legend deterministic and stable across samples); palette categories come
    first in their config order, any leftover columns follow sorted by name.
    """
    cols = list(ct.columns)
    if isinstance(color_map, dict) and color_map:
        in_palette = [c for c in color_map.keys() if c in cols]
        rest = [c for c in cols if c not in in_palette]
        return in_palette + sorted(rest, key=str)
    return cols


# ── core drawing ────────────────────────────────────────────────────────────

def draw_stacked_composition(ax, ct, color_map=None, *, normalize=False,
                             ylabel=None, xlabel=None, title=None,
                             legend_title=None, legend=True,
                             legend_kwargs=None, edge_lw=EDGE_LW):
    """Draw a stacked composition barplot onto an existing Axes.

    Parameters
    ----------
    ax : matplotlib Axes
    ct : DataFrame  (index = x groups, columns = stacked categories)
    color_map : dict {category -> hex} or None
    normalize : if True rows are rescaled to 100 (relative %)
    """
    color_map = color_map if isinstance(color_map, dict) else {}
    order = _ordered_columns(ct, color_map)          # order[0] = TOP of stack
    palette = _resolve_palette(order, color_map)     # config wins; rest distinct
    data = ct[order].copy()

    if normalize:
        row_sum = data.sum(axis=1).replace(0, pd.NA)
        data = data.div(row_sum, axis=0) * 100

    # Stacked bars draw the first column at the BOTTOM. Plot reversed so that
    # order[0] ends up on TOP, matching the legend built below.
    plot_cols = order[::-1]
    colors = [palette[str(c)] for c in plot_cols]
    data[plot_cols].plot(
        kind="bar", stacked=True, ax=ax, color=colors,
        edgecolor=EDGE_COLOR, linewidth=edge_lw, width=0.85, legend=False,
    )

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if title is not None:
        ax.set_title(title)
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(45)
        lbl.set_ha("right")

    if legend:
        # Legend top->bottom == order[0..] == stack top->bottom.
        handles = [Patch(facecolor=palette[str(c)],
                         edgecolor=EDGE_COLOR, linewidth=edge_lw, label=str(c))
                   for c in order]
        kw = dict(title=legend_title, bbox_to_anchor=(1.02, 1.0),
                  loc="upper left", fontsize=7, frameon=False)
        if legend_kwargs:
            kw.update(legend_kwargs)
        ax.legend(handles=handles, **kw)
    return ax


def save_stacked_composition(ct, out_path, color_map=None, *, normalize=False,
                             ylabel=None, xlabel=None, title=None,
                             legend_title=None, figsize=None, dpi=300,
                             edge_lw=EDGE_LW):
    """Standalone-PNG wrapper around :func:`draw_stacked_composition`."""
    if ct is None or len(ct) == 0:
        return
    if figsize is None:
        figsize = (max(8, len(ct) * 1.2), 6)
    fig, ax = plt.subplots(figsize=figsize)
    draw_stacked_composition(ax, ct, color_map, normalize=normalize,
                             ylabel=ylabel, xlabel=xlabel, title=title,
                             legend_title=legend_title, edge_lw=edge_lw)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def composition_pair(adata, group_key, cat_key, color_map, out_dir, prefix,
                     *, group_label=None, cat_label=None, dpi=300):
    """Write the absolute + relative PNG pair for crosstab(group_key x cat_key).

    `cat_key` columns are the stacked segments (use the annotation column so the
    stacks are cell types and colours come from annotation_colors[cat_key]).
    """
    os.makedirs(out_dir, exist_ok=True)
    if group_key not in adata.obs.columns or cat_key not in adata.obs.columns:
        return
    ct = pd.crosstab(adata.obs[group_key], adata.obs[cat_key])
    if ct.empty:
        return
    group_label = group_label or group_key
    cat_label = cat_label or cat_key.replace("cell_type_", "")
    save_stacked_composition(
        ct, os.path.join(out_dir, f"{prefix}_absolute.png"), color_map,
        normalize=False, ylabel="Number of cells", xlabel=group_label,
        title=f"{cat_label} by {group_label} (absolute)",
        legend_title=cat_label, dpi=dpi)
    save_stacked_composition(
        ct, os.path.join(out_dir, f"{prefix}_relative.png"), color_map,
        normalize=True, ylabel="Percentage (%)", xlabel=group_label,
        title=f"{cat_label} by {group_label} (relative)",
        legend_title=cat_label, dpi=dpi)


# ── sample × region composition (two multi-dimensional layouts) ──────────────

def _composition_long(obs, sample_key, region_key, cat_key, normalize):
    """Tidy counts/percent of cat_key within each (sample, region)."""
    df = (obs.groupby([sample_key, region_key, cat_key], observed=True)
             .size().reset_index(name="count"))
    if normalize:
        df["value"] = (df.groupby([sample_key, region_key], observed=True)["count"]
                       .transform(lambda x: 100 * x / x.sum() if x.sum() else x * 0))
    else:
        df["value"] = df["count"].astype(float)
    return df


def composition_grouped(adata, outer_key, inner_key, cat_key, color_map,
                        out_path, *, normalize=True, outer_order=None,
                        inner_order=None, cat_label="Cell type", ylabel=None,
                        dpi=300):
    """Single-axis stacked barplot: x = (outer × inner), grouped by `outer_key`
    with a bold header above each group + a dashed separator, and `inner_key`
    labels on the x-axis. Stacks = cat_key categories. Figure width scales with
    the number of columns so bars/labels stay readable."""
    obs = adata.obs
    if any(k not in obs.columns for k in (outer_key, inner_key, cat_key)):
        return
    df = _composition_long(obs, outer_key, inner_key, cat_key, normalize)
    if df.empty:
        return
    for key, ordering in ((outer_key, outer_order), (inner_key, inner_order)):
        if ordering:
            present = [v for v in ordering if v in set(df[key].astype(str))]
            if present:
                df[key] = pd.Categorical(df[key].astype(str),
                                         categories=present, ordered=True)
    pivot = (df.pivot_table(index=[outer_key, inner_key], columns=cat_key,
                            values="value", fill_value=0, observed=True)
               .sort_index(level=[0, 1]))
    if pivot.empty:
        return
    order = _ordered_columns(pivot, color_map)         # top->bottom
    palette = _resolve_palette(order, color_map)
    plot_cols = order[::-1]

    n_cols = len(pivot.index)
    fig_w = max(10, min(80, n_cols * 0.5 + 3))         # scale, but bounded
    fig, ax = plt.subplots(figsize=(fig_w, 6.5))
    pivot[plot_cols].plot(kind="bar", stacked=True, ax=ax,
                          color=[palette[str(c)] for c in plot_cols],
                          edgecolor=EDGE_COLOR, linewidth=0.3, width=0.85,
                          legend=False)
    # Inner labels vertical → never overlap regardless of count or name length.
    ax.set_xticklabels([str(inner) for _, inner in pivot.index],
                       rotation=90, ha="center", va="top", fontsize=6)
    ax.tick_params(axis="x", length=0)
    ax.margins(x=0.005)

    # Bold outer-group headers on two alternating vertical levels (so adjacent
    # sample names never collide) + dashed separators between groups.
    outer_positions = {}
    for i, (outer, _) in enumerate(pivot.index):
        outer_positions.setdefault(outer, []).append(i)
    top = 100.0 if normalize else ax.get_ylim()[1]
    hdr_fs = 8 if len(outer_positions) <= 12 else 7
    last = None
    for k, (outer, idxs) in enumerate(outer_positions.items()):
        y = top * (1.045 if (k % 2) else 1.01)          # stagger two levels
        ax.text((idxs[0] + idxs[-1]) / 2, y, str(outer),
                ha="center", va="bottom", fontsize=hdr_fs, fontweight="bold",
                clip_on=False)
        if last is not None:
            ax.axvline(idxs[0] - 0.5, color="gray", linestyle="--", linewidth=0.5)
        last = idxs[-1]

    ax.set_ylabel(ylabel or ("Percentage (%)" if normalize else "Number of cells"))
    ax.set_xlabel("")
    if normalize:
        ax.set_ylim(0, 100)
    handles = [Patch(facecolor=palette[str(c)], edgecolor=EDGE_COLOR,
                     linewidth=0.3, label=str(c)) for c in order]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              frameon=False, title=cat_label, ncol=1, fontsize=8, title_fontsize=9)
    fig.tight_layout(rect=[0, 0, 0.88, 1])
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def find_niche_column(adata, niche_column=None):
    """Return a usable spatial-niche obs column, or None if not available.

    Checks an explicit config name first, then a few conventional names.
    A column qualifies only if present with more than one distinct value.
    """
    candidates = []
    if niche_column:
        candidates.append(niche_column)
    candidates += ["niche", "spatial_niche", "niche_annotation", "cell_niche"]
    seen = set()
    for c in candidates:
        if not c or c in seen:
            continue
        seen.add(c)
        if c in adata.obs.columns:
            try:
                if adata.obs[c].nunique(dropna=True) > 1:
                    return c
            except Exception:
                continue
    return None
