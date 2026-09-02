#!/usr/bin/env python3
"""Rasterize the TACO workflow figure directly, at exact pixel geometry.

Drawn rather than converted: the figure is only rectangles, straight lines and
text, so drawing it gives exact control and avoids the aspect distortion an
SVG thumbnailer introduces. Same grid and palette as make_workflow_figure.py.

Usage:  python3 render_png.py <out.png> [scale]
"""
import sys
from PIL import Image, ImageDraw, ImageFont

REG = "/System/Library/Fonts/Supplemental/Arial.ttf"
BLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

S = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0   # supersample factor
GRID = 8 * S
def g(n): return int(round(n * GRID))

W, H = g(160), g(141)
LANE_X, LANE_W = g(15), g(84)
GATE_X, GATE_W = g(107), g(46)
BAND_X, BAND_W = g(4), g(5)
ROW_H, PITCH, SUB_H = g(6), g(9), g(5)

INK, MUTE = (0x14, 0x23, 0x2B), (0x5A, 0x6B, 0x75)
STYLE = {
    "data":    ((0xEF, 0xF2, 0xF4), (0x5A, 0x6B, 0x75)),
    "measure": ((0xFF, 0xFF, 0xFF), (0x00, 0x81, 0x5F)),
    "modify":  ((0xDC, 0xE9, 0xF2), (0x0B, 0x6F, 0xA4)),
    "gate":    ((0xFB, 0xE7, 0xD8), (0xC2, 0x54, 0x00)),
}
LW = max(1, int(round(1.6 * S)))

im = Image.new("RGB", (W, H), (255, 255, 255))
d = ImageDraw.Draw(im)
_fc = {}
def F(size, bold=False):
    k = (size, bold)
    if k not in _fc:
        _fc[k] = ImageFont.truetype(BLD if bold else REG, int(round(size * S)))
    return _fc[k]

def ctext(x, y, t, size, bold=False, col=INK, anchor="mm"):
    d.text((x, y), t, font=F(size, bold), fill=col, anchor=anchor)

def box(x, y, w, h, kind, lines, size=11.5, bold_first=True):
    fill, stroke = STYLE[kind]
    d.rectangle([x, y, x + w, y + h], fill=fill, outline=stroke, width=LW)
    lead = (size + 3.5) * S
    y0 = y + h / 2 - (len(lines) - 1) * lead / 2
    for i, t in enumerate(lines):
        ctext(x + w / 2, y0 + i * lead, t,
              size if i == 0 else size - 1.5,
              bold=(i == 0 and bold_first),
              col=INK if i == 0 else MUTE)

def vlink(x, y1, y2):
    hd = int(6 * S)
    d.line([x, y1, x, y2 - hd], fill=INK, width=LW)
    d.polygon([(x - 4 * S, y2 - hd), (x + 4 * S, y2 - hd), (x, y2)], fill=INK)

def hlink(x1, x2, y):
    hd = int(6 * S)
    d.line([x1, y, x2 - hd, y], fill=INK, width=LW)
    d.polygon([(x2 - hd, y - 4 * S), (x2 - hd, y + 4 * S), (x2, y)], fill=INK)

def band(y1, y2, label, sub):
    d.rectangle([BAND_X, y1, BAND_X + BAND_W, y2], fill=(255, 255, 255),
                outline=MUTE, width=max(1, int(1.4 * S)))
    cy = (y1 + y2) / 2
    for dx, txt, size, bold, col in ((-5 * S, label, 12.5, True, INK),
                                     (11 * S, sub, 10, False, MUTE)):
        t = Image.new("RGBA", (int(y2 - y1), int(g(3))), (0, 0, 0, 0))
        ImageDraw.Draw(t).text((t.width / 2, t.height / 2), txt,
                               font=F(size, bold), fill=col, anchor="mm")
        t = t.rotate(90, expand=True)
        im.paste(t, (int(BAND_X + BAND_W / 2 + dx - t.width / 2),
                     int(cy - t.height / 2)), t)

ctext(BAND_X, g(3), "TACO — run nine long-read assemblers, compare them on one table",
      15, bold=True, anchor="ls")

y = g(6); rows = []; cx = LANE_X + LANE_W / 2
a_top = y - g(1)
box(LANE_X, y, LANE_W, ROW_H, "data",
    ["One long-read set", "PacBio HiFi · Nanopore · PacBio CLR"]); rows.append((y, ROW_H)); y += PITCH
box(LANE_X, y, LANE_W, ROW_H, "measure",
    ["Install and select",
     "one conda environment; assemblers absent from PATH are skipped by name"]); rows.append((y, ROW_H)); y += PITCH

asms = ["HiCanu", "NextDenovo", "Peregrine", "IPA", "Flye", "hifiasm", "LJA", "MBG", "Raven"]
gap = g(1); bw = (LANE_W - gap * (len(asms) - 1)) // len(asms)
ctext(LANE_X, y - 6 * S, "Run each assembler on the same reads, in its own directory",
      10.5, col=MUTE, anchor="ls")
for i, a in enumerate(asms):
    box(LANE_X + i * (bw + gap), y, bw, SUB_H, "modify", [a], size=9.5)
rows.append((y, SUB_H)); y += PITCH - g(1)

box(LANE_X, y, LANE_W, ROW_H, "modify",
    ["Standardize",
     "uniform contig names; capture alternate contig sets the assemblers already wrote"]); rows.append((y, ROW_H)); y += PITCH

mets = ["BUSCO", "QUAST", "Telomere classification", "Merqury"]
mw = (LANE_W - gap * (len(mets) - 1)) // len(mets)
ctext(LANE_X, y - 6 * S, "Measure every assembly the same way", 10.5, col=MUTE, anchor="ls")
for i, m in enumerate(mets):
    box(LANE_X + i * (mw + gap), y, mw, SUB_H, "measure", [m], size=10)
rows.append((y, SUB_H)); y += PITCH - g(1)

box(LANE_X, y, LANE_W, ROW_H, "data",
    ["assembly_info.csv — the comparison table",
     "one row per metric, one column per assembly"]); rows.append((y, ROW_H))
band(a_top, y + ROW_H + g(2), "MEASURE", "the primary output")
y += PITCH + g(1)

ry = y - g(3)
d.line([BAND_X, ry, W - g(4), ry], fill=INK, width=max(2, int(2.4 * S)))
d.rectangle([LANE_X, ry - g(1.5), LANE_X + g(37), ry + g(1.5)], fill=(255, 255, 255))
ctext(LANE_X + g(1), ry, "--assembly-only stops here", 11.5, bold=True, anchor="lm")

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
    rows.append((y, ROW_H)); y += PITCH

box(LANE_X, y, LANE_W, ROW_H, "data",
    ["Final assembly + provenance + comparison report"]); rows.append((y, ROW_H))
band(b_top, y + ROW_H + g(2), "MODIFY", "optional refinement")

for i in range(len(rows) - 1):
    y1, h1 = rows[i]
    vlink(cx, y1 + h1, rows[i + 1][0])

ly = H - g(4); lx = LANE_X
for kind, label in (("measure", "measures only"), ("modify", "changes sequence"),
                    ("gate", "check that can veto the step beside it"),
                    ("data", "file on disk")):
    fill, stroke = STYLE[kind]
    d.rectangle([lx, ly - 9 * S, lx + g(2), ly + g(1.5) - 9 * S],
                fill=fill, outline=stroke, width=max(1, int(1.4 * S)))
    ctext(lx + g(2.6), ly, label, 10.5, col=MUTE, anchor="lm")
    lx += g(2.6) + int(len(label) * 5.6 * S) + g(3)

out = sys.argv[1] if len(sys.argv) > 1 else "taco_workflow.png"
im = im.resize((int(W / S * 2), int(H / S * 2)), Image.LANCZOS)   # 2x final
im.save(out, "PNG", optimize=True)
print(f"wrote {out}  {im.size[0]}x{im.size[1]}")
