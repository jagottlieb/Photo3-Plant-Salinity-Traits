"""Digitize the Na+ (Figure 1) and Cl- (Figure 2) bar plots in Campos-Villarreal et al. (2017).

The figures are embedded grayscale raster images, so values are read from pixels:
gray y-axis ticks calibrate each panel, bar tops give the means, and the top cap of
each black error bar gives mean + SE.

Usage:
    python digitize_campos_villarreal.py               # extract, write CSV + QC overlays to output/
    python digitize_campos_villarreal.py --write-xlsx  # also write the sheet into the workbook
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
from PIL import Image, ImageDraw

from workbook_io import write_sheets

LIT = Path(__file__).resolve().parents[1]
PDF = LIT / "Campos-Villareal 2017 (1).pdf"
SHEET = "Campos-Villareal Ion Accum"
OUT = LIT.parents[1] / "output" / "campos_villarreal_digitized"

CULTIVARS = ["Riverside", "SVA1", "Apache", "SVA2"]  # bar order within each salinity group
OUTLINED = {"Riverside", "Apache"}                   # bars drawn with a black outline
SALINITY = [0.08, 2.5, 3.0, 3.5]                     # dS/m, group order along x
COMPARTMENTS = ["leaf", "stem", "root"]              # panels a, b, c from top to bottom

FIGURES = [
    # page is 0-based; image is the index in page.get_images(); ticks are the y-axis labels per panel
    dict(figure=1, ion="Na+", page=4, image=0,
         ticks=[np.arange(0, 1601, 200), np.arange(0, 1601, 200), np.arange(0, 1001, 200)]),
    dict(figure=2, ion="Cl-", page=5, image=0,
         ticks=[np.arange(0, 81, 10), np.arange(0, 61, 10), np.arange(0, 61, 10)]),
]
UNITS = "mg/kg dry weight"

INK = 200    # gray axes/ticks and patterned fills are darker than this
BLACK = 90   # bar outlines and error bars are darker than this


def runs(mask):
    """(start, stop) index pairs of consecutive True values."""
    m = np.concatenate([[0], np.asarray(mask, int), [0]])
    d = np.diff(m)
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def load_image(fig):
    doc = pymupdf.open(PDF)
    xref = doc[fig["page"]].get_images(full=True)[fig["image"]][0]
    pix = pymupdf.Pixmap(doc, xref)
    return np.array(Image.frombytes("L", (pix.width, pix.height), pix.samples)).astype(int)


def find_panels(a):
    """(ax0, ax1, top, bottom) per panel: y-axis column span and its row span."""
    ink = a < INK
    xaxes = [(s, e) for s, e in runs(ink.mean(axis=1) > 0.5)]
    panels, prev = [], 0
    for s, e in xaxes:
        region = ink[prev:e]
        frac = region.mean(axis=0)
        peak = int(np.argmax(frac[: a.shape[1] // 2]))
        ax0, ax1 = next((c0, c1 - 1) for c0, c1 in runs(frac > 0.8 * frac[peak]) if c0 <= peak < c1)
        top, bottom = max(runs(region[:, (ax0 + ax1) // 2]), key=lambda r: r[1] - r[0])
        panels.append((ax0, ax1, prev + top, prev + bottom))
        prev = e
    return panels


def find_ticks(a, ax0, top, bottom):
    """Row centers of the tick marks just left of the y-axis."""
    band = a[top - 5:bottom + 5, ax0 - 8:ax0 - 2] < INK
    rows = band.all(axis=1)
    return np.array([top - 5 + (s + e - 1) / 2 for s, e in runs(rows)])


def find_bars(a, ax1, baseline):
    """x-extents of the 16 bars, from the band of rows just above the x-axis."""
    band = a[baseline - 40:baseline - 8, ax1 + 3:] < INK
    cols = band.mean(axis=0) > 0.3
    for s, e in runs(~cols):  # bridge 1-2 px gaps in the dotted SVA2 fill
        if e - s <= 2 and s > 0:
            cols[s:e] = True
    r = [(s + ax1 + 3, e + ax1 + 3) for s, e in runs(cols)]
    widths = np.array([e - s for s, e in r])
    w = np.median(widths[widths > 15])
    bars, i = [], 0
    while i < len(r) and len(bars) < 16:
        s, e = r[i]
        if e - s > w / 2:
            bars.append((s, e))
            i += 1
            continue
        # white Apache bars show up as thin runs (outline, error-bar stem, outline) spanning one bar width
        j = i
        while j + 1 < len(r) and r[j + 1][1] - s <= 1.25 * w:
            j += 1
        if abs(r[j][1] - s - w) < w / 4:
            bars.append((s, r[j][1]))
        i = j + 1
    return bars, w


def outline_width(a, bars, baseline):
    """Thickness of the black bar outlines, from the left edge of the Apache bars."""
    ws = []
    for s, e in bars[2::4]:
        row = a[baseline - 20, s:e] < BLACK
        ws.append(runs(row)[0][1] - runs(row)[0][0])
    return float(np.median(ws))


def bar_top(a, s, e, baseline, outlined, t):
    """Pixel row of the bar's value: top edge of the fill, or center of the top outline."""
    w = e - s
    cols = np.r_[s + int(0.12 * w):s + int(0.28 * w), e - int(0.28 * w):e - int(0.12 * w)]
    ink = (a[:, cols] < INK).mean(axis=1)
    r = baseline - 3
    while r > 0 and ink[r - 6:r].max() > 0.05:
        r -= 1
    top = r
    while ink[top] < 0.05:
        top += 1
    return top + (t - 1) / 2 if outlined else top - 0.5


def error_cap(a, s, e, top):
    """Pixel row of the center of the error-bar's upper cap (None if hidden)."""
    xc = (s + e) // 2
    win = a[:, xc - 6:xc + 7] < BLACK
    stem = xc - 6 + int(np.argmax(win[int(top) - 15:int(top) - 2].sum(axis=0)))
    r = int(top) - 2
    gap = 0
    last = None
    while r > 0 and gap < 3:
        if (a[r, stem - 1:stem + 2] < BLACK).any():
            last, gap = r, 0
        else:
            gap += 1
        r -= 1
    if last is None:
        return None
    w = e - s
    capw = (a[last - 2:last + 12, int(s + 0.25 * w):int(e - 0.25 * w)] < BLACK).mean(axis=1)
    cap_rows = np.where(capw > 0.6)[0] + last - 2
    if len(cap_rows) == 0:
        return last
    first = runs(capw > 0.6)[0]
    return last - 2 + (first[0] + first[1] - 1) / 2


def digitize_figure(fig):
    a = load_image(fig)
    panels = find_panels(a)
    assert len(panels) == 3, f"Figure {fig['figure']}: found {len(panels)} panels"
    rows, overlay = [], Image.fromarray(a.astype(np.uint8)).convert("RGB")
    draw = ImageDraw.Draw(overlay)
    for p, ((ax0, ax1, top, bottom), labels, comp) in enumerate(zip(panels, fig["ticks"], COMPARTMENTS)):
        ticks = find_ticks(a, ax0, top, bottom)
        assert len(ticks) == len(labels), f"Fig {fig['figure']}{'abc'[p]}: {len(ticks)} ticks vs {len(labels)} labels"
        slope, icpt = np.polyfit(ticks, labels[::-1], 1)
        resid = np.abs(np.polyval([slope, icpt], ticks) - labels[::-1]).max()
        to_val = lambda y: slope * y + icpt
        baseline = int(round(ticks[-1]))
        bars, w = find_bars(a, ax1, baseline)
        assert len(bars) == 16, f"Fig {fig['figure']}{'abc'[p]}: found {len(bars)} bars"
        t = outline_width(a, bars, baseline)
        print(f"Fig {fig['figure']}{'abc'[p]}: {abs(slope):.3f} units/px, tick fit resid {resid:.2f}, "
              f"bar width {w:.0f}px, outline {t:.0f}px")
        for i, (s, e) in enumerate(bars):
            cv, sal = CULTIVARS[i % 4], SALINITY[i // 4]
            y_mean = bar_top(a, s, e, baseline, cv in OUTLINED, t)
            y_cap = error_cap(a, s, e, y_mean)
            mean = to_val(y_mean)
            se = to_val(y_cap) - mean if y_cap is not None else np.nan
            rows.append(dict(cultivar=cv, salinity_dS_m=sal, compartment=comp, ion=fig["ion"],
                             value=round(mean, 2), se=round(se, 2), units=UNITS,
                             error_type="SE (n=240 obs)", figure=fig["figure"], panel="abc"[p],
                             page=fig["page"] + 1, method="raster pixel extraction"))
            draw.line([(s, y_mean), (e, y_mean)], fill=(255, 0, 0), width=2)
            if y_cap is not None:
                xc = (s + e) / 2
                draw.line([(xc - w / 3, y_cap), (xc + w / 3, y_cap)], fill=(0, 160, 255), width=2)
        for y in ticks:
            draw.line([(ax0 - 20, y), (ax0 - 2, y)], fill=(0, 200, 0), width=2)
    OUT.mkdir(parents=True, exist_ok=True)
    overlay.save(OUT / f"qc_figure{fig['figure']}.png")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-xlsx", action="store_true")
    args = ap.parse_args()
    df = pd.DataFrame([r for fig in FIGURES for r in digitize_figure(fig)])
    df.to_csv(OUT / "campos_villarreal_ions.csv", index=False)
    print(df.pivot_table(index=["ion", "compartment", "salinity_dS_m"], columns="cultivar",
                         values="value", sort=False).round(1).to_string())
    if args.write_xlsx:
        notes = [f"Digitized from Campos-Villarreal et al. 2017 Figures 1 (Na+) and 2 (Cl-) "
                 f"with data/literature_data/salt_accumulation_analysis/{Path(__file__).name}",
                 "Means over n=240 observations; error bars are SE of the mean."]
        write_sheets([(SHEET, notes, df)])


if __name__ == "__main__":
    main()
