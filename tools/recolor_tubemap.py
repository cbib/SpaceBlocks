#!/usr/bin/env python3
"""Recolour a Snakevision tube-map SVG by SpaceBlocks block and add a legend.

Uses a high-contrast, colour-blind-friendly categorical palette. Head rules are
identified by their suffix; core rules are assigned explicitly to a phase.
Unmapped rules are coloured grey and reported.

Usage:
    python recolor_tubemap.py rulegraph.svg [output.svg]
"""
from pathlib import Path
import re
import sys

HEAD = {
    "_vhd": ("#000000", "Visium HD"),
    "_x5k": ("#E69F00", "Xenium 5K"),
    "_ate": ("#56B4E9", "Atera"),
    "_mer": ("#009E73", "MERSCOPE"),
}

CORE = {
    "pre": ("#D55E00", "Preprocessing"),
    "post": ("#0072B2", "Postprocessing"),
    "explore": ("#CC79A7", "Exploration"),
}

UNKNOWN = "#777777"

CORE_PHASE = {
    # Preprocessing
    "validate_input": "pre",
    "qc_sweep": "pre",
    "qc_sweep_all": "pre",
    "preprocess_umap": "pre",
    "leiden_analysis": "pre",
    "generate_annotation_template": "pre",
    "ingest_ref": "pre",
    "spatial_niches": "pre",
    "run_preprocessing": "pre",
    # Postprocessing
    "annotate_cells": "post",
    "integrate_samples": "post",
    "pseudobulk_aggregate": "post",
    "pseudobulk_de": "post",
    "neighbourhood_analysis": "post",
    "subcluster": "post",
    "subcluster_all": "post",
    "sample_report": "post",
    "run_postprocessing": "post",
    # Exploration
    "explore_genes": "explore",
    "explore_genes_integrated": "explore",
    "explore_genes_sample": "explore",
    "run_exploration": "explore",
}

def colour_for(rule):
    """Return the block colour for a rule, or None when it is unmapped."""
    if rule is None:
        return None
    for suffix, (colour, _) in HEAD.items():
        if rule.endswith(suffix):
            return colour
    phase = CORE_PHASE.get(rule)
    return CORE[phase][0] if phase else None


def recolour(svg):
    nodes = {}
    for match in re.finditer(
        r'<circle cx="([\d.]+)" cy="([\d.]+)" fill="#[0-9a-fA-F]{6}" id="N(\w+)"',
        svg,
    ):
        nodes[(float(match.group(1)), float(match.group(2)))] = match.group(3)

    segments = []
    for match in re.finditer(
        r'<line stroke="#[0-9a-fA-F]{6}" stroke-width="2.0" '
        r'x1="([\d.]+)" x2="([\d.]+)" y1="([\d.]+)" y2="([\d.]+)"',
        svg,
    ):
        x1, x2, y1, y2 = map(float, match.groups())
        if x1 == x2 and (x1, y1) in nodes:
            segments.append((x1, min(y1, y2), max(y1, y2), nodes[(x1, y1)]))

    def owner_at(x, y):
        for segment_x, low, high, rule in segments:
            if abs(segment_x - x) < 0.01 and low - 0.01 <= y <= high + 0.01:
                return rule
        return None

    output = []
    last_path = None
    unknown = set()

    for token in re.finditer(r'<[^>]+>|[^<]+', svg):
        element = token.group(0)
        node_match = re.search(r'id="N(\w+)"', element)

        if element.startswith("<circle") and node_match:
            rule = node_match.group(1)
            colour = colour_for(rule)
            if colour is None:
                unknown.add(rule)
                colour = UNKNOWN
            element = re.sub(r'fill="#[0-9a-fA-F]{6}"', f'fill="{colour}"', element)

        elif element.startswith("<path") and 'stroke="#' in element:
            move_match = re.search(r'M ([\d.]+),([\d.]+)', element)
            last_path = (
                owner_at(float(move_match.group(1)), float(move_match.group(2)))
                if move_match else None
            )
            colour = colour_for(last_path)
            if colour:
                element = re.sub(
                    r'stroke="#[0-9a-fA-F]{6}"', f'stroke="{colour}"', element
                )

        elif element.startswith("<line") and 'stroke="#' in element:
            coords = re.search(
                r'x1="([\d.]+)" x2="([\d.]+)" y1="([\d.]+)" y2="([\d.]+)"',
                element,
            )
            if coords:
                x1, x2, y1, _ = map(float, coords.groups())
                rule = (
                    nodes.get((x1, y1)) or owner_at(x1, y1)
                    if x1 == x2 else last_path
                )
                colour = colour_for(rule)
                if colour:
                    element = re.sub(
                        r'stroke="#[0-9a-fA-F]{6}"', f'stroke="{colour}"', element
                    )

        output.append(element)

    if unknown:
        sys.stderr.write(
            "[warn] unmapped rules coloured grey: "
            + ", ".join(sorted(unknown))
            + "\n"
        )
    return "".join(output)


def add_legend(svg):
    match = re.search(
        r'viewBox="([\d.]+)[,\s]([\d.]+)[,\s]([\d.]+)[,\s]([\d.]+)"', svg
    )
    if not match:
        sys.stderr.write("[warn] no viewBox found; legend skipped\n")
        return svg

    min_x, min_y, width, height = map(float, match.groups())
    row_height = 22
    heading_height = 24
    padding = 24
    legend_height = padding + heading_height + max(len(HEAD), len(CORE)) * row_height + 12
    legend_y = min_y + height + padding

    svg = svg.replace(
        match.group(0),
        f'viewBox="{min_x},{min_y},{width},{height + legend_height}"',
    )

    left_x = min_x + 10
    right_x = min_x + max(220, width * 0.52)

    def entry(x, y, colour, label):
        return (
            f'<circle cx="{x}" cy="{y}" r="6" fill="{colour}" '
            f'stroke="white" stroke-width="2"/>'
            f'<text x="{x + 13}" y="{y + 4}" font-family="sans-serif" '
            f'font-size="14">{label}</text>'
        )

    parts = [
        '<g id="legend">',
        f'<line x1="{min_x + 10}" x2="{min_x + width - 10}" '
        f'y1="{legend_y - 10}" y2="{legend_y - 10}" stroke="#D0D0D0"/>',
        f'<text x="{left_x}" y="{legend_y + 6}" font-family="sans-serif" '
        f'font-size="15" font-weight="bold">Headblocks</text>',
        f'<text x="{right_x}" y="{legend_y + 6}" font-family="sans-serif" '
        f'font-size="15" font-weight="bold">Coreblocks</text>',
    ]

    for index, (_, (colour, label)) in enumerate(HEAD.items()):
        parts.append(entry(left_x + 7, legend_y + 28 + index * row_height, colour, label))

    for index, (_, (colour, label)) in enumerate(CORE.items()):
        parts.append(entry(right_x + 7, legend_y + 28 + index * row_height, colour, label))

    parts.append("</g>")
    return svg.replace("</svg>", "".join(parts) + "</svg>")


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: recolor_tubemap.py input.svg [output.svg]")

    input_path = Path(sys.argv[1])
    output_path = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else input_path.with_name(f"{input_path.stem}_recoloured.svg")
    )

    svg = input_path.read_text(encoding="utf-8")
    output_path.write_text(add_legend(recolour(svg)), encoding="utf-8")
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
