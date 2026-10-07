"""Digitize the Na+ (Figure 1) and K+ (Figure 2) bar plots in Sanchez-Ledesma et al. (2022).

Unlike Campos-Villarreal, these figures are vector graphics. The bars are pattern fills that
PyMuPDF does not list, but every bar has a vector error bar: two vertical segments that meet at
the bar top (the mean) and end in caps at mean +/- SE. The y-axis tick marks calibrate the values
and the x-axis category ticks assign each error bar to a treatment.

Usage:
    python digitize_sanchez_ledesma.py               # extract, write CSV + QC overlays to output/
    python digitize_sanchez_ledesma.py --write-xlsx  # also write the sheet into the workbook
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
from PIL import Image, ImageDraw

from workbook_io import write_sheets

LIT = Path(__file__).resolve().parents[1]
PDF = next(LIT.glob("S*Ledesma*2022*.pdf"))
SHEET = "Sanchez-Ledesma Ion Accum"
OUT = LIT.parents[1] / "output" / "sanchez_ledesma_digitized"

# (label, NaCl mM, inoculated with Scleroderma sp.) in x order
TREATMENTS = [("Testigo", 0, False), ("0 mM", 0, True), ("20 mM", 20, True),
              ("25 mM", 25, True), ("30 mM", 30, True), ("35 mM", 35, True)]
COMPARTMENTS = ["root", "stem", "leaf"]  # bar order within each treatment (legend: Raiz, Tallo, Hoja)

FIGURES = [
    # page is 0-based; ticks are the y-axis labels from bottom to top
    dict(figure=1, ion="Na+", page=5, ticks=np.arange(0, 71, 10)),
    dict(figure=2, ion="K+", page=6, ticks=np.arange(0, 46, 5)),
]
UNITS = "mg/g dry weight"
QC_DPI = 300


def line_items(d):
    """((x0, y0), (x1, y1)) for each straight segment of a drawing."""
    return [((it[1].x, it[1].y), (it[2].x, it[2].y)) for it in d["items"] if it[0] == "l"]


def find_axes(drawings):
    """y-axis tick rows (bottom to top) and x-axis category boundaries (left to right)."""
    yticks = xticks = None
    for d in drawings:
        segs = line_items(d)
        if d["type"] != "s" or len(segs) < 3 or len(segs) != len(d["items"]):
            continue
        if all(abs(a[1] - b[1]) < 0.01 and abs(a[0] - b[0]) < 4 for a, b in segs):
            yticks = sorted((a[1] for a, b in segs), reverse=True)
        elif all(abs(a[0] - b[0]) < 0.01 and abs(a[1] - b[1]) < 4 for a, b in segs):
            xticks = sorted(a[0] for a, b in segs)
    return np.array(yticks), np.array(xticks)


def find_error_bars(drawings):
    """(x, y_mean, y_upper_cap, y_lower_cap) for each error bar."""
    bars = []
    for d in drawings:
        segs = line_items(d)
        if d["type"] != "fs" or len(segs) != 4:
            continue
        vert = [s for s in segs if abs(s[0][0] - s[1][0]) < 0.01]
        if len(vert) != 2 or vert[0][0] != vert[1][0]:
            continue
        x, y_mean = vert[0][0]
        ends = sorted(v[1][1] for v in vert)
        bars.append((x, y_mean, ends[0], ends[1]))
    return sorted(bars)


def digitize_figure(doc, fig):
    page = doc[fig["page"]]
    drawings = page.get_drawings()
    yticks, xticks = find_axes(drawings)
    labels = fig["ticks"]
    assert len(yticks) == len(labels), f"Fig {fig['figure']}: {len(yticks)} ticks vs {len(labels)} labels"
    assert len(xticks) == len(TREATMENTS) + 1, f"Fig {fig['figure']}: {len(xticks)} x-axis ticks"
    slope, icpt = np.polyfit(yticks, labels, 1)
    resid = np.abs(np.polyval([slope, icpt], yticks) - labels).max()
    to_val = lambda y: slope * y + icpt
    bars = find_error_bars(drawings)
    assert len(bars) == len(TREATMENTS) * len(COMPARTMENTS), f"Fig {fig['figure']}: {len(bars)} error bars"
    print(f"Fig {fig['figure']} ({fig['ion']}): {abs(slope):.4f} units/pt, tick fit resid {resid:.4f}")

    scale = QC_DPI / 72
    pix = page.get_pixmap(dpi=QC_DPI)
    overlay = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    draw = ImageDraw.Draw(overlay)
    rows = []
    for x, y_mean, y_up, y_lo in bars:
        t = int(np.searchsorted(xticks, x)) - 1
        in_group = sorted(b[0] for b in bars if xticks[t] < b[0] < xticks[t + 1])
        comp = COMPARTMENTS[in_group.index(x)]
        label, nacl, inoc = TREATMENTS[t]
        mean = to_val(y_mean)
        assert abs((y_mean - y_up) - (y_lo - y_mean)) < 0.15, f"asymmetric error bar at x={x}"
        se = (to_val(y_up) - to_val(y_lo)) / 2
        rows.append(dict(treatment=label, nacl_mM=nacl, inoculated=inoc, compartment=comp,
                         ion=fig["ion"], value=round(mean, 3), se=round(se, 3), units=UNITS,
                         error_type="SE (n=5)", figure=fig["figure"], page=fig["page"] + 1,
                         method="vector error-bar extraction"))
        draw.line([((x - 4) * scale, y_mean * scale), ((x + 4) * scale, y_mean * scale)], fill=(255, 0, 0), width=3)
        draw.line([((x - 3) * scale, y_up * scale), ((x + 3) * scale, y_up * scale)], fill=(0, 160, 255), width=3)
    x0, x1 = xticks[0] - 45, xticks[-1] + 15
    y0, y1 = yticks[-1] - 15, yticks[0] + 40
    OUT.mkdir(parents=True, exist_ok=True)
    overlay.crop(tuple(int(v * scale) for v in (x0, y0, x1, y1))).save(OUT / f"qc_figure{fig['figure']}.png")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-xlsx", action="store_true")
    args = ap.parse_args()
    doc = pymupdf.open(PDF)
    df = pd.DataFrame([r for fig in FIGURES for r in digitize_figure(doc, fig)])
    df.to_csv(OUT / "sanchez_ledesma_ions.csv", index=False)
    wide = df.pivot_table(index=["ion", "treatment"], columns="compartment", values="value", sort=False)
    print(wide[COMPARTMENTS].round(2).to_string())
    if args.write_xlsx:
        notes = [f"Digitized from Sanchez-Ledesma et al. 2022 Figures 1 (Na+) and 2 (K+) "
                 f"with data/literature_data/salt_accumulation_analysis/{Path(__file__).name}",
                 "Testigo = not inoculated, no NaCl; all other treatments inoculated with Scleroderma sp."]
        write_sheets([(SHEET, notes, df)])


if __name__ == "__main__":
    main()
