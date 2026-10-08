"""First day after burn-in, med run: s = 0.45 (batch 20261002_1436) vs s = 0.32 (batch 20261002_1521).

Both batches used 1 m compartments and ran before the gp lag fix. Shows why E is lower at
s = 0.32: soil-root conductance comparable to plant conductance lowers psi_b, stem storage
and pre-dawn psi_l, so stomata start partly closed.
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

repo = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
sys.path.insert(0, str(repo / "sensitivity_analysis"))
import post_processing as pp

res = repo / "sensitivity_analysis/Phase 1 - pi0/results"
runs = {"s = 0.45": res / "20261002_1436/timeseries/med__S2.csv",
        "s = 0.32": res / "20261002_1521/timeseries/med__S2.csv"}
LAI, F_CAP = 17, 0.5


def pair_sum(col):
    return np.array([sum(json.loads(v)) for v in col])


data = {}
for name, path in runs.items():
    d = pd.read_csv(path)
    d = d[(d.days_since_burn_in > 0) & (d.days_since_burn_in < 1 - 1e-9)].copy()
    d["hour"] = (d.days_since_burn_in % 1) * 24
    d["gsr_tot"] = pair_sum(d["gsr"])
    d["qs_tot"] = pair_sum(d["qs"])
    d["gpl"] = d["gp"] * LAI
    d["g_series"] = 1 / (1 / d["gsr_tot"] + 1 / d["gpl"])
    data[name] = d

cols = ["ev", "qs_tot", "qw_stem", "qw_leaf", "w_stem", "psi_w_stem", "psi_l", "psi_x", "psi_b",
        "gsw", "gpl", "gsr_tot", "g_series"]
summary = pd.DataFrame({
    name: {**{f"max {c}": d[c].max() for c in cols}, **{f"min {c}": d[c].min() for c in cols},
           "daily E (mm)": d["ev"].sum() * 1800 / 1000, "daily root uptake (mm)": d["qs_tot"].sum() * 1800 / 1000,
           "daily net stem storage (mm)": d["qw_stem"].sum() * 1800 / 1000}
    for name, d in data.items()})
pd.set_option("display.width", 200)
print(summary.round(4).to_string())

# Midday snapshot at peak E
for name, d in data.items():
    r = d.loc[d["ev"].idxmax()]
    print(f"\n{name} at peak E (hour {r.hour:.1f}): E={r.ev:.4f} qs={r.qs_tot:.4f} qw_stem={r.qw_stem:.4f} "
          f"psi_s~{json.loads(r.psi_s)[0]:.3f} psi_b={r.psi_b:.3f} psi_x={r.psi_x:.3f} psi_l={r.psi_l:.3f} "
          f"gsw={r.gsw:.4f} gp*LAI={r.gpl:.4f} gsr={r.gsr_tot:.3f} w_stem={r.w_stem:.3f} psi_w_stem={r.psi_w_stem:.3f}")

pp.set_plot_style(10)
panels = [("ev", "E"), ("qs_tot", "Root uptake $q_s$"), ("qw_stem", "Stem storage flux"),
          ("w_stem", "Stem storage fraction"), ("psi_l", "$\\psi_l$"), ("psi_b", "$\\psi_b$ (root base)"),
          ("gsw", "$g_{sw}$"), ("gpl", "$g_p\\cdot$LAI")]
fig, axes = plt.subplots(4, 2, figsize=(11, 10), sharex=True)
for ax, (c, lab) in zip(axes.ravel(), panels):
    for name, d in data.items():
        ax.plot(d["hour"], d[c], label=name)
    ax.set_title(lab, fontsize="medium")
axes[0, 0].legend(frameon=False)
for ax in axes[-1]:
    ax.set_xlabel("Hour of day (first day after burn-in)")
    ax.set_xticks(range(0, 25, 6))
fig.suptitle("med pi0: s = 0.45 vs s = 0.32, first day after burn-in")
fig.tight_layout()
fig.savefig(Path(__file__).resolve().parent / "compare_s045_s032_day1.png", dpi=130, bbox_inches="tight")
print("saved")
