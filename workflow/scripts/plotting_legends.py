"""Helpers for keeping categorical legends inside multi-panel figures."""

from __future__ import annotations

import math
import textwrap

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_LEGEND_FONTSIZE = 12


def category_labels(adata, obs_key):
    """Return displayed category labels in their plotting order."""
    if not obs_key or obs_key not in adata.obs.columns:
        return []
    values = adata.obs[obs_key]
    categories = getattr(values.dtype, "categories", None)
    if categories is not None:
        labels = categories
    else:
        labels = values.dropna().unique()
    return [str(label) for label in labels]


def _wrap_label(label, width=24):
    """Wrap only the displayed legend text; underlying category names are unchanged."""
    return textwrap.fill(
        str(label),
        width=width,
        break_long_words=True,
        break_on_hyphens=False,
    )


def legend_layout(labels, wrap_width=24, preferred_rows=5, max_columns=4,
                  fontsize=DEFAULT_LEGEND_FONTSIZE):
    """Return wrapped labels, column count, and estimated height in inches."""
    wrapped = [_wrap_label(label, wrap_width) for label in labels]
    if not wrapped:
        return wrapped, 1, 0.0

    ncol = min(max_columns, max(1, math.ceil(len(wrapped) / preferred_rows)))
    rows = math.ceil(len(wrapped) / ncol)
    line_counts = [label.count("\n") + 1 for label in wrapped]

    # Matplotlib fills legends column-first. Estimate each displayed row using the
    # tallest wrapped entry that can land on it, then reserve a title/margin band.
    row_lines = []
    for row in range(rows):
        row_lines.append(
            max(
                (line_counts[col * rows + row]
                 for col in range(ncol)
                 if col * rows + row < len(line_counts)),
                default=1,
            )
        )
    # Reserve physical space in proportion to the rendered font size. The
    # multiplier includes Matplotlib's default line spacing plus row padding.
    line_height = 0.30 * (fontsize / 10)
    title_height = 0.55 * (fontsize / 10)
    height = title_height + sum(line_height * lines for lines in row_lines)
    return wrapped, ncol, max(0.8, height)


def legend_height(labels_by_panel, minimum=0.8):
    """Height needed by the largest panel legend."""
    heights = [legend_layout(labels)[2] for labels in labels_by_panel]
    return max([minimum, *heights])


def panel_figure(n_panels, labels_by_panel, plot_rows=1, panel_width=8.0,
                 plot_height=6.0, wspace=0.18, hspace=0.25,
                 top_margin=0.0):
    """Create plot axes plus a dedicated legend row in the same figure.

    The returned legend axes are part of the figure itself. Consequently, a PDF
    writer always keeps each legend on the same page as its plots.
    """
    legend_h = legend_height(labels_by_panel)
    fig = plt.figure(
        figsize=(panel_width * n_panels,
                 plot_height * plot_rows + legend_h + top_margin)
    )
    grid = fig.add_gridspec(
        plot_rows + 1,
        n_panels,
        height_ratios=[plot_height] * plot_rows + [legend_h],
        wspace=wspace,
        hspace=hspace,
    )
    plot_axes = np.asarray(
        [[fig.add_subplot(grid[row, col]) for col in range(n_panels)]
         for row in range(plot_rows)]
    )
    legend_axes = [fig.add_subplot(grid[plot_rows, col]) for col in range(n_panels)]
    for ax in legend_axes:
        ax.set_axis_off()
    return fig, plot_axes, legend_axes


def grid_figure(nrows, ncols, legend_labels, cell_width=5.0, cell_height=4.0,
                wspace=0.18, hspace=0.25, top_margin=0.0,
                max_legend_columns=8):
    """Create a plot grid with one shared legend band on the same figure."""
    _, _, legend_h = legend_layout(
        legend_labels,
        max_columns=max_legend_columns,
    )
    legend_h = max(0.8, legend_h)
    fig = plt.figure(
        figsize=(cell_width * ncols,
                 cell_height * nrows + legend_h + top_margin)
    )
    grid = fig.add_gridspec(
        nrows + 1,
        ncols,
        height_ratios=[cell_height] * nrows + [legend_h],
        wspace=wspace,
        hspace=hspace,
    )
    axes = np.asarray(
        [[fig.add_subplot(grid[row, col]) for col in range(ncols)]
         for row in range(nrows)]
    )
    legend_ax = fig.add_subplot(grid[nrows, :])
    legend_ax.set_axis_off()
    return fig, axes, legend_ax


def draw_legend_on_axis(
    legend_ax,
    handles,
    labels,
    title=None,
    fontsize=DEFAULT_LEGEND_FONTSIZE,
    wrap_width=24,
    preferred_rows=5,
    max_columns=4,
):
    """Draw handles in a reserved legend axis using adaptive rows and wrapping."""
    legend_ax.set_axis_off()
    if not handles:
        return
    wrapped, ncol, _ = legend_layout(
        labels,
        wrap_width=wrap_width,
        preferred_rows=preferred_rows,
        max_columns=max_columns,
        fontsize=fontsize,
    )
    legend_ax.legend(
        handles,
        wrapped,
        title=title,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=ncol,
        fontsize=fontsize,
        title_fontsize=fontsize + 1,
        frameon=False,
        handletextpad=0.35,
        columnspacing=0.8,
        labelspacing=0.35,
        borderaxespad=0.0,
    )


def move_legend_to_axis(
    source_axes,
    legend_ax,
    title=None,
    fontsize=DEFAULT_LEGEND_FONTSIZE,
    wrap_width=24,
    preferred_rows=5,
    max_columns=4,
):
    """Move unique entries from Scanpy legends into a dedicated in-figure axis.

    Any duplicate legends on the supplied source axes are removed. ``source_axes``
    can therefore contain matching UMAP and spatial panels while only one legend
    is rendered below them.
    """
    if not isinstance(source_axes, (list, tuple, np.ndarray)):
        source_axes = [source_axes]

    handles = []
    labels = []
    seen_labels = set()
    for ax in source_axes:
        legend = ax.get_legend()
        if legend is None:
            continue
        axis_handles = (
            list(legend.legend_handles)
            if hasattr(legend, "legend_handles")
            else list(getattr(legend, "legendHandles", []))
        )
        axis_labels = [text.get_text() for text in legend.get_texts()]
        for handle, label in zip(axis_handles, axis_labels):
            if label not in seen_labels:
                handles.append(handle)
                labels.append(label)
                seen_labels.add(label)
        legend.remove()

    draw_legend_on_axis(
        legend_ax,
        handles,
        labels,
        title=title,
        fontsize=fontsize,
        wrap_width=wrap_width,
        preferred_rows=preferred_rows,
        max_columns=max_columns,
    )
