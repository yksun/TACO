#!/usr/bin/env python3
"""Render the TACO workflow diagram as SVG on a strict grid.

Deliberate constraints, so the output reads as a drawn vector figure rather than
a generated one:
  * rectangles, straight lines and straight arrows only -- no rounded cards,
    circles, curves, gradients, shadows or 3-D;
  * every coordinate is a multiple of GRID;
  * boxes at the same level share exact dimensions;
  * connectors run strictly vertical or strictly horizontal.

Usage:  python3 make_workflow_figure.py <out.svg>
"""
import sys

GRID = 8
def g(n): return n * GRID                      # grid units -> px

# ── canvas ──────────────────────────────────────────────────────────────────
W, H = g(160), g(141)
LANE_X, LANE_W = g(15), g(84)                  # main flow column
GATE_X, GATE_W = g(107), g(46)                 # gate column
BAND_X, BAND_W = g(4), g(5)                    # band label bar
ROW_H, PITCH = g(6), g(9)                      # box height, row pitch
SUB_H = g(5)

INK, MUTE, RULE = "#14232B", "#5A6B75", "#14232B"
STYLE = {                                       # fill, stroke
    "data":    ("#EFF2F4", "#5A6B75"),
    "measure": ("#FFFFFF", "#00815F"),
    "modify":  ("#DCE9F2", "#0B6FA4"),
    "gate":    ("#FBE7D8", "#C25400"),
}
o = []


def esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def box(x, y, w, h, kind, lines, size=11.5, bold_first=True):
    fill, stroke = STYLE[kind]
    o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" '
             f'stroke="{stroke}" stroke-width="1.6"/>')
    n = len(lines)
    lead = size + 3.5
    y0 = y + h / 2 - (n - 1) * lead / 2 + size * 0.36
    for i, t in enumerate(lines):
        weight = "600" if (i == 0 and bold_first) else "400"
        fs = size if i == 0 else size - 1.5
        col = INK if i == 0 else MUTE
        o.append(f'<text x="{x + w/2}" y="{y0 + i*lead:.0f}" font-size="{fs}" '
                 f'font-weight="{weight}" fill="{col}" text-anchor="middle">'
                 f'{esc(t)}</text>')


def vlink(x, y1, y2, head=True):
    o.append(f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2-(6 if head else 0)}" '
             f'stroke="{INK}" stroke-width="1.6"/>')
    if head:
        o.append(f'<polygon points="{x-4},{y2-6} {x+4},{y2-6} {x},{y2}" fill="{INK}"/>')


def hlink(x1, x2, y, head=True):
    o.append(f'<line x1="{x1}" y1="{y}" x2="{x2-(6 if head else 0)}" y2="{y}" '
             f'stroke="{INK}" stroke-width="1.6"/>')
    if head:
        o.append(f'<polygon points="{x2-6},{y-4} {x2-6},{y+4} {x2},{y}" fill="{INK}"/>')


def band(y1, y2, label, sub):
    o.append(f'<rect x="{BAND_X}" y="{y1}" width="{BAND_W}" height="{y2-y1}" '
             f'fill="#FFFFFF" stroke="{MUTE}" stroke-width="1.4"/>')
    cy = (y1 + y2) / 2
    o.append(f'<text x="{BAND_X + BAND_W/2 - 5}" y="{cy}" font-size="12.5" '
             f'font-weight="700" fill="{INK}" text-anchor="middle" '
             f'transform="rotate(-90 {BAND_X + BAND_W/2 - 5} {cy})">{esc(label)}</text>')
    o.append(f'<text x="{BAND_X + BAND_W/2 + 11}" y="{cy}" font-size="10" '
             f'fill="{MUTE}" text-anchor="middle" '
             f'transform="rotate(-90 {BAND_X + BAND_W/2 + 11} {cy})">{esc(sub)}</text>')


def main(out):
    cx = LANE_X + LANE_W / 2
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
             f'viewBox="0 0 {W} {H}" font-family="Helvetica Neue,Helvetica,Arial,sans-serif">')
    o.append(f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>')
    o.append(f'<text x="{BAND_X}" y="{g(3)}" font-size="15" font-weight="700" '
             f'fill="{INK}">TACO — run nine long-read assemblers, compare them on one table</text>')

    y = g(6)
    rows = []

    # ── band A: measure ─────────────────────────────────────────────────────
    a_top = y - g(1)
    box(LANE_X, y, LANE_W, ROW_H, "data",
        ["One long-read set", "PacBio HiFi · Nanopore · PacBio CLR"]); rows.append(y); y += PITCH
    box(LANE_X, y, LANE_W, ROW_H, "measure",
        ["Install and select", "one conda environment; assemblers absent from PATH are skipped by name"]); rows.append(y); y += PITCH

    asms = ["HiCanu", "NextDenovo", "Peregrine", "IPA", "Flye",
            "hifiasm", "LJA", "MBG", "Raven"]
    n = len(asms)
    gap = g(1)
    bw = (LANE_W - gap * (n - 1)) // n
    o.append(f'<text x="{LANE_X}" y="{y-6}" font-size="10.5" fill="{MUTE}">'
             f'Run each assembler on the same reads, in its own directory</text>')
    for i, a in enumerate(asms):
        box(LANE_X + i * (bw + gap), y, bw, SUB_H, "modify", [a], size=9.5)
    rows.append(y); y += PITCH - g(1)

    box(LANE_X, y, LANE_W, ROW_H, "modify",
        ["Standardize", "uniform contig names; capture alternate contig sets the assemblers already wrote"]); rows.append(y); y += PITCH

    mets = ["BUSCO", "QUAST", "Telomere classification", "Merqury"]
    gap2 = g(1)
    mw = (LANE_W - gap2 * (len(mets) - 1)) // len(mets)
    o.append(f'<text x="{LANE_X}" y="{y-6}" font-size="10.5" fill="{MUTE}">'
             f'Measure every assembly the same way</text>')
    for i, m in enumerate(mets):
        box(LANE_X + i * (mw + gap2), y, mw, SUB_H, "measure", [m], size=10)
    rows.append(y); y += PITCH - g(1)

    box(LANE_X, y, LANE_W, ROW_H, "data",
        ["assembly_info.csv — the comparison table",
         "one row per metric, one column per assembly"]); rows.append(y)
    a_bot = y + ROW_H + g(2)
    band(a_top, a_bot, "MEASURE", "the primary output")
    y += PITCH + g(1)

    # ── the dividing rule ───────────────────────────────────────────────────
    ry = y - g(3)
    o.append(f'<line x1="{BAND_X}" y1="{ry}" x2="{W-g(4)}" y2="{ry}" '
             f'stroke="{RULE}" stroke-width="2.4"/>')
    o.append(f'<rect x="{LANE_X}" y="{ry-g(1.5)}" width="{g(37)}" height="{g(3)}" fill="#FFFFFF"/>')
    o.append(f'<text x="{LANE_X+g(1)}" y="{ry+4}" font-size="11.5" font-weight="700" '
             f'fill="{INK}">--assembly-only stops here</text>')

    # ── band B: modify ──────────────────────────────────────────────────────
    b_top = y - g(1)
    steps = [
        ("measure", ["Select a backbone", "explicit score, written to disk with its weights"], None),
        ("modify",  ["Build the telomere pool", "telomere-bearing contigs from every assembly"], None),
        ("modify",  ["Merge", "rebuild from the pool, dropping covered backbone sequence"],
                    ["GATE", "restore backbone contigs", "carrying genes the merge dropped"]),
        ("modify",  ["purge_dups", "remove residual duplication"],
                    ["GATE", "accept only within the", "taxon BUSCO tolerance"]),
        ("modify",  ["Screen foreign sequence", "read depth AND composition both required"], None),
        ("modify",  ["Resolve mis-joins", "split where other assemblies disagree and no read spans the junction"],
                    ["GATE", "majority of informative assemblies", "must split, and reads must not span"]),
        ("measure", ["Measure what was delivered", "on the file TACO actually writes"],
                    ["GATE", "report backbone vs delivered:", "size, telomeres, gene content"]),
    ]
    for kind, lines, gate in steps:
        box(LANE_X, y, LANE_W, ROW_H, kind, lines)
        if gate:
            box(GATE_X, y, GATE_W, ROW_H, "gate", gate, size=10.5)
            hlink(LANE_X + LANE_W, GATE_X, y + ROW_H / 2)
        rows.append(y); y += PITCH

    box(LANE_X, y, LANE_W, ROW_H, "data",
        ["Final assembly + provenance + comparison report"])
    rows.append(y)
    b_bot = y + ROW_H + g(2)
    band(b_top, b_bot, "MODIFY", "optional refinement")

    # ── vertical connectors, drawn last so they sit under nothing ───────────
    for i in range(len(rows) - 1):
        y1 = rows[i] + (SUB_H if rows[i] in () else ROW_H)
        h1 = SUB_H if i in (2, 4) else ROW_H
        vlink(cx, rows[i] + h1, rows[i + 1])

    # ── legend, aligned to the same grid ────────────────────────────────────
    ly = H - g(4)
    lx = LANE_X
    for kind, label in (("measure", "measures only"), ("modify", "changes sequence"),
                        ("gate", "check that can veto the step beside it"),
                        ("data", "file on disk")):
        fill, stroke = STYLE[kind]
        o.append(f'<rect x="{lx}" y="{ly-9}" width="{g(2)}" height="{g(1.5)}" '
                 f'fill="{fill}" stroke="{stroke}" stroke-width="1.4"/>')
        o.append(f'<text x="{lx+g(2.6)}" y="{ly+1}" font-size="10.5" fill="{MUTE}">{esc(label)}</text>')
        lx += g(2.6) + len(label) * 5.6 + g(3)
    o.append("</svg>")
    open(out, "w").write("\n".join(o))
    print(f"wrote {out}  ({W}x{H})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "workflow.svg")
