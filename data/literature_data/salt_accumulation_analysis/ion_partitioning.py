"""Treatment-averaged ion concentrations per compartment, with compartment ratios, for each study.

Reads the digitized '<study> Ion Accum' sheets (see digitize_*.py) and writes a '<study> Ion Ratios'
sheet per study. Ratios are given on absolute concentrations and, for the salt ions, on the increase
above the no-salt control, which isolates where the added salt went.

With --latex, the salt-ion concentrations are also converted to mM in tissue water and written as
two LaTeX tables (concentrations; compartment ratios) at the assumed water contents +/- WC_DELTA.

Usage:
    python ion_partitioning.py               # print tables
    python ion_partitioning.py --write-xlsx  # also write the sheets into the workbook
    python ion_partitioning.py --latex       # also write ion_tables.tex next to this script
"""
import argparse
from pathlib import Path

import pandas as pd

from workbook_io import XLSX, write_sheets

COMPARTMENTS = ["leaf", "stem", "root"]
RATIOS = [("leaf", "stem"), ("leaf", "root"), ("stem", "root")]
SALT_IONS = {"Na+", "Cl-"}
MW = {"Na+": 22.99, "Cl-": 35.45, "K+": 39.10}          # g/mol
MG_PER_G = {"mg/kg dry weight": 1e-3, "mg/g dry weight": 1.0}
# Fresh-mass water content (water / fresh mass). Leaf: lab pecan leaf samples (Leaf Data tab).
# Stem: saturated sapwood from storage_volume_calcs.py, 0.61 m3 water per 600 kg dry wood.
# Root: no measurements; assumed slightly wetter than stem.
WATER_CONTENT = {"leaf": 0.61, "stem": 0.50, "root": 0.60}
WC_DELTA = 0.10                                          # +/- absolute change in water content
TEX = Path(__file__).resolve().parent / "ion_tables.tex"
# LaTeX tables: Campos-Villareal Na+ only. Its Cl- axis is labelled mg/kg (implausibly low; units
# uncertain); Sanchez-Ledesma reports ~4% leaf Na+ even in controls, unrealistic for pecan.
TEX_STUDY, TEX_ION = "Campos-Villareal", "Na+"

STUDIES = [
    # level: treatment column; control: baseline for 'above control';
    # no_salt: levels left out of the salt-treatment mean and the 'above control' rows
    # tex_level: LaTeX label format for a treatment level
    dict(name="Campos-Villareal", level="salinity_dS_m", control=0.08, no_salt=[0.08],
         tex_level="{} dS m$^{{-1}}$", note="Means over the 4 cultivars."),
    dict(name="Sanchez-Ledesma", level="treatment", control="0 mM", no_salt=["Testigo", "0 mM"],
         tex_level="{}",
         note="Testigo = not inoculated, no NaCl; other treatments inoculated with Scleroderma sp.; "
              "'above control' uses the inoculated 0 mM treatment."),
]


def compartment_means(ions: pd.DataFrame, level: str) -> pd.DataFrame:
    """Mean over replicates/cultivars, one row per ion and level (in sheet order), one column per compartment."""
    ions = ions.copy()
    for col in ["ion", level]:
        ions[col] = pd.Categorical(ions[col], categories=list(dict.fromkeys(ions[col])))
    wide = ions.pivot_table(index=["ion", level], columns="compartment", values="value",
                            aggfunc="mean", observed=True)
    return wide[COMPARTMENTS]


def partition_table(ions: pd.DataFrame, study: dict) -> pd.DataFrame:
    """Absolute and above-control rows per ion, each followed by the mean over the salt treatments."""
    level, control = study["level"], study["control"]
    means = compartment_means(ions, level)
    rows = []
    for basis in ["absolute", "above control"]:
        for ion in means.index.get_level_values("ion").unique():
            if basis == "above control" and ion not in SALT_IONS:
                continue
            sub = means.loc[ion]
            salt = [lv for lv in sub.index if lv not in study["no_salt"]]
            shown = list(sub.index) if basis == "absolute" else salt
            vals = sub - sub.loc[control] if basis == "above control" else sub
            for lv in shown:
                rows.append({"ion": ion, level: str(lv), "basis": basis, **vals.loc[lv]})
            rows.append({"ion": ion, level: f"{salt[0]}-{salt[-1]} mean", "basis": basis,
                         **vals.loc[salt].mean()})
    return add_ratios(pd.DataFrame(rows))


def add_ratios(table: pd.DataFrame) -> pd.DataFrame:
    for num, den in RATIOS:
        table[f"{num}/{den}"] = table[num] / table[den]
    table["leaf share"] = table["leaf"] / table[COMPARTMENTS].sum(axis=1)
    return table


def to_millimolar(table: pd.DataFrame, units: str, delta: float = 0.0) -> pd.DataFrame:
    """Dry-mass concentrations -> mM in tissue water at WATER_CONTENT + delta; ratios recomputed."""
    out = table.copy()
    for c in COMPARTMENTS:
        wc = WATER_CONTENT[c] + delta
        out[c] = table[c] * MG_PER_G[units] * 1000 / table["ion"].map(MW) / (wc / (1 - wc))
    return add_ratios(out)


def tex_tabular(header: list[list[str]], rows: list[list[str]], colspec: str) -> list[str]:
    """booktabs tabular with a \\midrule under each header row; a row of None draws a \\midrule."""
    line = lambda r: r"\midrule" if r is None else " & ".join(r) + r" \\"
    return [rf"\begin{{tabular}}{{{colspec}}}", r"\toprule", *(f"{line(h)}\n\\midrule" for h in header),
            *map(line, rows), r"\bottomrule", r"\end{tabular}"]


def latex_tables(study: dict, table: pd.DataFrame, units: str) -> str:
    """One table: TEX_ION concentrations (range at +/- WC_DELTA) above leaf/stem, leaf/root ratios."""
    def rng(mid, a, b):
        f = ".0f" if abs(mid) >= 10 else ".1f" if abs(mid) >= 1 else ".2f"
        return f"{mid:{f}} ({min(a, b):{f}}--{max(a, b):{f}})"
    def label(lv):
        return "Salt-treatment mean" if "mean" in lv else study["tex_level"].format(lv)
    level = study["level"]
    mid, wet, dry = (to_millimolar(table, units, d) for d in (0.0, WC_DELTA, -WC_DELTA))
    ion = table["ion"] == TEX_ION
    conc, ratio = [], []
    for i in table.index[ion & (table["basis"] == "absolute")]:
        lv = table.at[i, level]
        conc += [None] * ("mean" in lv) + [[label(lv)] + [rng(mid.at[i, c], wet.at[i, c], dry.at[i, c])
                                                          for c in COMPARTMENTS]]
        above = table.index[ion & (table["basis"] == "above control") & (table[level] == lv)]
        vals = [f"{mid.at[k, r]:.2f}" for k in [i, *above] for r in ("leaf/stem", "leaf/root")]
        ratio += [None] * ("mean" in lv) + [[label(lv)] + (vals + ["--", "--"])[:4]]
    wc = ", ".join(f"{c} {100 * w:.0f}\\%" for c, w in WATER_CONTENT.items())
    caption = (r"Na$^+$ in pecan rootstocks from \citet{CamposVillarreal2017} (means over four cultivars), "
               f"converted to concentrations in tissue water at fresh-mass water contents of {wc}. Top: "
               f"concentrations (mM); parentheses give the range for water contents $\\pm${100 * WC_DELTA:.0f} "
               "percentage points. Bottom: leaf-to-stem (L/S) and leaf-to-root (L/R) concentration ratios, on "
               "absolute concentrations and on the increase above the 0.08~dS~m$^{-1}$ control.")
    conc_tab = tex_tabular([["Concentration (mM)", "Leaf", "Stem", "Root"]], conc, "lrrr")
    ratio_tab = tex_tabular([["Concentration ratios", r"\multicolumn{2}{c}{Absolute}",
                              r"\multicolumn{2}{c}{Above control}"],
                             ["Treatment", "L/S", "L/R", "L/S", "L/R"]], ratio, "lrrrr")
    return "\n".join([r"\begin{table}[htbp]", r"\centering", r"\small", rf"\caption{{{caption}}}",
                      r"\label{tab:na_partitioning}", *conc_tab, r"\par\medskip", *ratio_tab,
                      r"\end{table}"]) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-xlsx", action="store_true")
    ap.add_argument("--latex", action="store_true")
    args = ap.parse_args()
    sheets, results = [], []
    for study in STUDIES:
        ions = pd.read_excel(XLSX, sheet_name=f"{study['name']} Ion Accum", skiprows=2)
        units = ions["units"].iloc[0]
        table = partition_table(ions, study)
        print(f"\n=== {study['name']} ({units})")
        with pd.option_context("display.width", 200, "display.max_columns", 20):
            print(table.round(2).to_string(index=False))
        notes = [f"Mean concentrations ({units}) and compartment ratios from '{study['name']} Ion Accum'; "
                 f"made with data/literature_data/salt_accumulation_analysis/{Path(__file__).name}",
                 study["note"] + " 'above control' = treatment minus control; "
                 "'leaf share' = leaf / (leaf + stem + root)"]
        sheets.append((f"{study['name']} Ion Ratios", notes, table.round(3)))
        results.append((study, table, units))
    if args.write_xlsx:
        write_sheets(sheets)
    if args.latex:
        TEX.parent.mkdir(parents=True, exist_ok=True)
        TEX.write_text(latex_tables(*next(r for r in results if r[0]["name"] == TEX_STUDY)), encoding="utf-8")
        print(f"\n{TEX.read_text(encoding='utf-8')}\nWrote {TEX}")


if __name__ == "__main__":
    main()
