#!/usr/bin/env python3
"""Recolour a snakevision tube-map SVG by SpaceBlocks block (viridis) and add a legend.

Snakevision auto-assigns node colours (recycling them across unrelated rules) and ignores
the Graphviz `color` attribute, so colouring has to be done on the rendered SVG. Each
element is mapped to its rule by POSITION — circles by their `id`, tube segments/curves by
the node they descend from — never by colour, since the same snakevision colour can land on
rules in different blocks.

Head rules are detected by their `_vhd` / `_x5k` suffix (so a new head just works); core
rules are mapped to a phase below. Unmapped rules are coloured grey and reported.

Usage:
    python recolor_tubemap.py input.svg [output.svg]
"""
import re
import sys

# ── Block palette (viridis) ──────────────────────────────────────────────────
HEAD = {  # rule-name suffix -> (hex, legend label)
    "_vhd": ("#440154", "Visium HD"),
    "_x5k": ("#414487", "Xenium 5K"),
}
CORE_PHASE = {  # core rule -> phase (names don't encode the phase, so map it here)
    "validate_input": "pre", "qc_sweep": "pre", "preprocess_umap": "pre",
    "leiden_analysis": "pre", "generate_annotation_template": "pre",
    "ingest_ref": "pre", "spatial_niches": "pre",
    "annotate_cells": "post", "integrate_samples": "post", "pseudobulk_aggregate": "post",
    "pseudobulk_de": "post", "neighbourhood_analysis": "post", "subcluster": "post",
    "sample_report": "post",
    "explore_genes_integrated": "explore", "explore_genes_sample": "explore",
}
CORE = {  # phase -> (hex, legend label)
    "pre": ("#2a788e", "Preprocessing"),
    "post": ("#35b779", "Postprocessing"),
    "explore": ("#7ad151", "Exploration"),
}
UNKNOWN = "#999999"


def colour_for(rule):
    """Block colour for a rule, or None if unmapped."""
    if rule is None:
        return None
    for suf, (col, _) in HEAD.items():
        if rule.endswith(suf):
            return col
    phase = CORE_PHASE.get(rule)
    return CORE[phase][0] if phase else None


def recolour(svg):
    # circles -> node positions
    nodes = {}
    for m in re.finditer(r'<circle cx="([\d.]+)" cy="([\d.]+)" fill="#[0-9a-fA-F]{6}" id="N(\w+)"', svg):
        nodes[(float(m.group(1)), float(m.group(2)))] = m.group(3)
    # vertical coloured lines -> tube segments (owner = node at the top of the segment)
    segs = []
    for m in re.finditer(r'<line stroke="#[0-9a-fA-F]{6}" stroke-width="2.0" '
                         r'x1="([\d.]+)" x2="([\d.]+)" y1="([\d.]+)" y2="([\d.]+)"', svg):
        x1, x2, y1, y2 = map(float, m.groups())
        if x1 == x2 and (x1, y1) in nodes:
            segs.append((x1, min(y1, y2), max(y1, y2), nodes[(x1, y1)]))

    def owner_at(x, y):
        for sx, lo, hi, rule in segs:
            if abs(sx - x) < 0.01 and lo - 0.01 <= y <= hi + 0.01:
                return rule
        return None

    out, last_path, unknown = [], None, set()
    for tok in re.finditer(r'<[^>]+>|[^<]+', svg):
        s = tok.group(0)
        mc = re.search(r'id="N(\w+)"', s)
        if s.startswith("<circle") and mc:                       # node
            c = colour_for(mc.group(1))
            if c is None:
                unknown.add(mc.group(1)); c = UNKNOWN
            s = re.sub(r'fill="#[0-9a-fA-F]{6}"', f'fill="{c}"', s)
        elif s.startswith("<path") and 'stroke="#' in s:          # tube corner
            mm = re.search(r'M ([\d.]+),([\d.]+)', s)
            last_path = owner_at(float(mm.group(1)), float(mm.group(2))) if mm else None
            c = colour_for(last_path)
            if c:
                s = re.sub(r'stroke="#[0-9a-fA-F]{6}"', f'stroke="{c}"', s)
        elif s.startswith("<line") and 'stroke="#' in s:          # tube segment / run
            mx = re.search(r'x1="([\d.]+)" x2="([\d.]+)" y1="([\d.]+)" y2="([\d.]+)"', s)
            x1, x2, y1, y2 = map(float, mx.groups())
            rule = (nodes.get((x1, y1)) or owner_at(x1, y1)) if x1 == x2 else last_path
            c = colour_for(rule)
            if c:
                s = re.sub(r'stroke="#[0-9a-fA-F]{6}"', f'stroke="{c}"', s)
        out.append(s)
    if unknown:
        sys.stderr.write("[warn] unmapped rules coloured grey: %s\n" % ", ".join(sorted(unknown)))
    return "".join(out)


def add_legend(svg):
    m = re.search(r'viewBox="([\d.]+)[,\s]([\d.]+)[,\s]([\d.]+)[,\s]([\d.]+)"', svg)
    if not m:
        sys.stderr.write("[warn] no viewBox found; legend skipped\n")
        return svg
    minx, miny, w, h = map(float, m.groups())
    y0 = miny + h + 20
    svg = svg.replace(m.group(0), f'viewBox="{minx},{miny},{w},{h + 120}"')
    xa, xb = minx + 10, minx + 205

    def entry(cx, cy, col, txt):
        return (f'<circle cx="{cx}" cy="{cy}" r="6" fill="{col}" stroke="white" stroke-width="2"/>'
                f'<text x="{cx + 13}" y="{cy + 4}" font-family="sans-serif" font-size="14">{txt}</text>')

    parts = [
        '<g id="legend">',
        f'<line x1="{minx + 10}" x2="{minx + w - 10}" y1="{y0 - 8}" y2="{y0 - 8}" stroke="lightgrey"/>',
        f'<text x="{xa}" y="{y0 + 8}" font-family="sans-serif" font-size="15" font-weight="bold">Headblocks</text>',
        entry(xa + 7, y0 + 28, HEAD["_vhd"][0], HEAD["_vhd"][1]),
        entry(xa + 7, y0 + 48, HEAD["_x5k"][0], HEAD["_x5k"][1]),
        f'<text x="{xb}" y="{y0 + 8}" font-family="sans-serif" font-size="15" font-weight="bold">Coreblocks</text>',
        entry(xb + 7, y0 + 28, CORE["pre"][0], CORE["pre"][1]),
        entry(xb + 7, y0 + 48, CORE["post"][0], CORE["post"][1]),
        entry(xb + 7, y0 + 68, CORE["explore"][0], CORE["explore"][1]),
        '</g>',
    ]
    return svg.replace("</svg>", "".join(parts) + "</svg>")


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: recolor_tubemap.py input.svg [output.svg]")
    inp = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else inp.rsplit(".svg", 1)[0] + "_recoloured.svg"
    svg = add_legend(recolour(open(inp).read()))
    open(out, "w").write(svg)
    print("wrote", out)


if __name__ == "__main__":
    main()
