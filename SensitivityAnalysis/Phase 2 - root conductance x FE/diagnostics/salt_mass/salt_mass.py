"""Salt mass sanity check for a Phase 2 batch: how much salt is transported, not just how concentrations change.

For each run in one scenario, over the analysis window (after burn-in):
  - water taken up by roots (positive qs only) and water lost to soil (negative qs, hydraulic redistribution)
  - salt arriving at the root surface: sum over compartments of max(qs, 0) * cs
  - salt admitted to the plant (model MW_uptake_total) and the implied filtration efficiency
    1 - admitted / arriving, plus the admitted sap concentration (admitted / water taken up)
  - split of admitted salt into stem and leaf storage, and the leaf share
  - admitted salt as a fraction of the soil salt stock (the soil's salt mass is held fixed by the model,
    so uptake never depletes it and excluded salt never accumulates)
  - c_leaf if all admitted salt had gone to leaf storage, for comparison with the 24.3 mol/m3 threshold
  - per-plant totals in mg NaCl, using ground area per plant = LA / LAI

Writes salt_mass_<batch>_<scenario>.csv and a figure of cumulative salt mass per run group next to this script.
Conditions: Phase 2 config (wet constant soil, 100 mM), batch given on the command line.

Usage (from the repo root):
    python "SensitivityAnalysis/Phase 2 - root conductance x FE/diagnostics/salt_mass/salt_mass.py" <batch_id> [scenario]
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
PHASE = HERE.parents[1]
REPO = PHASE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "SensitivityAnalysis"))

import species_traits
from post_processing import run_label, set_plot_style
from soil import Berger

NACL_G_PER_MOL = 58.44

batch_id = sys.argv[1]
scenario = sys.argv[2] if len(sys.argv) > 2 else "S5"
batch = PHASE / "results" / batch_id
cfg = yaml.safe_load((batch / "config.yaml").read_text(encoding="utf-8"))
common = cfg["plots"]["common"]
species = getattr(species_traits, cfg["simulation"]["species"])()
lai, la = float(species.LAI), float(species.LA)
plant_ground_m2 = la / lai
porosity = Berger.N
zr = float(species.ZR)  # zr_arr: species
dt = float(cfg["simulation"]["timestepM"]) * 60.0

manifest = pd.read_csv(batch / "manifest.csv")
manifest = manifest[manifest["scenario"] == scenario]


def pairs(df: pd.DataFrame, col: str) -> np.ndarray:
    """Per-compartment column as an (n_steps, n_comp) array."""
    return np.array([json.loads(c) for c in df[col]])


rows, curves = [], {}
for r in manifest.itertuples():
    df = pd.read_csv(batch / r.timeseries_file)
    w = df[~df["is_burn_in"]]
    qs, cs, s = pairs(w, "qs"), pairs(w, "cs"), pairs(w, "s")
    q_in = np.clip(qs, 0, None)

    arriving_step = (q_in * cs).sum(axis=1) * 1e-6 * dt              # mol/m2 ground per step
    arriving = arriving_step.cumsum()
    base = {k: float(df[k].iloc[len(df) - len(w) - 1]) for k in ("MW_uptake_total", "MW_uptake_stem", "MW_uptake_leaf")}
    admitted = w["MW_uptake_total"].to_numpy() - base["MW_uptake_total"]
    to_stem = w["MW_uptake_stem"].to_numpy() - base["MW_uptake_stem"]
    to_leaf = w["MW_uptake_leaf"].to_numpy() - base["MW_uptake_leaf"]
    curves[r.run] = (w["days_since_burn_in"].to_numpy(), arriving, admitted, to_leaf, r)

    water_in_mm = q_in.sum() * dt * 1e-3                               # um/s * s -> um -> mm
    water_out_mm = -np.clip(qs, None, 0).sum() * dt * 1e-3
    soil_salt = float((cs[0] * zr * porosity * s[0]).sum())           # mol/m2 ground, all compartments
    leaf_vol = float(w["vw_leaf"].iloc[-1]) * lai                      # m3 water per m2 ground
    rows.append({
        "run": r.run,
        "group": r.group,
        "E_set": json.loads(r.overrides).get("E", cfg["baseline"]["E"]),
        "water_in_mm": water_in_mm,
        "water_out_mm": water_out_mm,
        "transpired_mm": w["ev"].sum() * dt * 1e-3,
        "salt_arriving_mol_m2": arriving[-1],
        "salt_admitted_mol_m2": admitted[-1],
        "implied_E": 1 - admitted[-1] / arriving[-1] if arriving[-1] > 0 else np.nan,
        "sap_conc_mol_m3": admitted[-1] / (water_in_mm * 1e-3) if water_in_mm > 0 else np.nan,
        "to_stem_mol_m2": to_stem[-1],
        "to_leaf_mol_m2": to_leaf[-1],
        "leaf_share": to_leaf[-1] / admitted[-1] if admitted[-1] > 0 else np.nan,
        "admitted_frac_of_soil_salt": admitted[-1] / soil_salt,
        "c_leaf_end": float(w["c_leaf"].iloc[-1]),
        "c_leaf_if_all_to_leaf": float(w["c_leaf"].iloc[0]) + admitted[-1] / leaf_vol,
        "admitted_mg_NaCl_per_plant": admitted[-1] * plant_ground_m2 * NACL_G_PER_MOL * 1e3,
        "to_leaf_mg_NaCl_per_plant": to_leaf[-1] * plant_ground_m2 * NACL_G_PER_MOL * 1e3,
    })

table = pd.DataFrame(rows)
table.to_csv(HERE / f"salt_mass_{batch_id}_{scenario}.csv", index=False)
pd.set_option("display.width", 200)
print(f"Batch {batch_id}, {scenario}; LAI {lai}, LA {la} m2 -> {plant_ground_m2 * 1e4:.2f} cm2 ground per plant;"
      f" soil salt stock {soil_salt:.2f} mol/m2")
print(table.drop(columns="group").to_string(index=False, float_format=lambda v: f"{v:.4g}"))

# Cumulative salt mass per run group: arriving at roots, admitted to plant, allocated to leaf.
set_plot_style(common.get("font_size", 10))
panels = [("Salt arriving at roots (mol m$^{-2}$)", 1), ("Salt admitted to plant (mol m$^{-2}$)", 2),
          ("Salt allocated to leaf (mmol m$^{-2}$)", 3)]
for group in table["group"].unique():
    fig, axes = plt.subplots(len(panels), 1, figsize=(8, 2.8 * len(panels)), sharex=True)
    for run, (x, *ys, r) in curves.items():
        if r.group != group:
            continue
        for ax, (label, k) in zip(axes, panels):
            y = ys[k - 1] * (1e3 if k == 3 else 1)
            ax.plot(x, y, color=common["run_colors"].get(run), ls=common["run_linestyles"].get(run, "-"),
                    lw=1.2, label=run_label(r.overrides, common.get("display_names", {})))
            ax.set_ylabel(label)
    axes[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize="small", frameon=False)
    axes[-1].set_xlabel("Days")
    fig.suptitle(f"Cumulative salt mass ({scenario}, {group})")
    fig.tight_layout()
    out = HERE / f"salt_mass_{batch_id}_{scenario}_{group}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.name}")
