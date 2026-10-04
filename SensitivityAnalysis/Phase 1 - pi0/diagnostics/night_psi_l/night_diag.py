"""Night psi_l diagnostic: why pre-dawn psi_l stays far below psi_s with 0.2 m compartments.

Uses batch 20261002_1550 (s_init 0.32, 0.2 m compartments, before the gp lag fix).
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
sys.path.insert(0, str(REPO / "SensitivityAnalysis"))
from post_processing import set_plot_style  # noqa: E402

CSV = REPO / "SensitivityAnalysis" / "Phase 1 - pi0" / "results" / "20261002_1550" / "timeseries" / "med__S2.csv"
HERE = Path(__file__).resolve().parent
LAI, F, VWT = 17.0, 0.5, 0.011

df = pd.read_csv(CSV)
pair = lambda c: np.array([json.loads(v) for v in df[c]])
df["psi_s1"] = pair("psi_s")[:, 0]
df["qs_sum"] = pair("qs").sum(1)
df["gsr_tot"] = pair("gsr").sum(1)
df["gpL2"] = 2 * df["gp"] * LAI          # each half of the xylem path: gp LAI / F and gp LAI / (1-F)
df["gwL"] = df["gw_stem"] * LAI
df["night"] = df["forcing_phi"] <= 0
# night n = the dark period that ends on the morning of analysis day n
t = df["days_since_burn_in"] + 0.25
df["night_id"] = np.floor(t - 1e-9).astype(int)
nights = df[df["night"] & (df["night_id"] >= 0) & ((t % 1) < 0.5)]
counts = nights.groupby("night_id").size()
nights = nights[nights["night_id"].isin(counts.index[counts == counts.max()])]
nm = nights.groupby("night_id").mean(numeric_only=True)
nm["tot_gradient"] = nm["psi_s1"] - nm["psi_l"]
nm["soil_drop"] = nm["psi_s1"] - nm["psi_b"]
nm["xylem_drop"] = nm["psi_b"] - nm["psi_l"]
nm["E_over_gsr"] = nm["ev"] / nm["gsr_tot"]

# stem storage capacitance (ground basis, um/MPa) from pre-dawn w_stem vs psi_w_stem
pd5 = df[np.isclose((df["time_days"] % 1) * 24, 5.0) & ~df["is_burn_in"]]
slope = np.polyfit(pd5["psi_w_stem"], pd5["w_stem"], 1)[0]
C = VWT * LAI * 1e6 * slope
nm["tau_days"] = C / nm["gsr_tot"] / 86400
night_len = df[df["night"]].shape[0] / (df.shape[0] / 48) * 1800  # s of darkness per day
nm["max_recharge_um"] = nm["gsr_tot"] * nm["soil_drop"] * night_len
nm["night_E_um"] = nm["ev"] * night_len

cols = ["psi_s1", "psi_b", "psi_x", "psi_w_stem", "psi_l", "soil_drop", "xylem_drop", "w_stem",
        "ev", "qs_sum", "qw_stem", "qw_leaf", "gsr_tot", "gpL2", "gwL", "tau_days", "max_recharge_um", "night_E_um"]
pd.set_option("display.width", 300, "display.precision", 4, "display.max_columns", 30)
print(f"stem capacitance C = {C:.0f} um/MPa (ground), dw/dpsi = {slope:.3f} /MPa; dark period = {night_len/3600:.1f} h")
print(nm.loc[[0, 1, 3, 5, 8, 11, 14], cols].to_string())

set_plot_style(10)
fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.5), sharex=True)
x = nm.index + 1
ax = axes[0]
for c, lab, col, ls in [("psi_s1", r"Soil $\psi_s$", "#8c6d31", "-"), ("psi_b", r"Root base $\psi_b$", "#bf812d", "--"),
                        ("psi_x", r"Stem node $\psi_x$", "#31a354", "--"), ("psi_w_stem", r"Stem storage $\psi_{w,stem}$", "#006d2c", ":"),
                        ("psi_l", r"Leaf $\psi_l$", "#08519c", "-")]:
    ax.plot(x, nm[c], ls, color=col, marker="o", ms=3, label=lab)
ax.set_ylabel("Night-mean potential (MPa)")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize="small")
ax.set_title("Night (dark-period mean), med run, batch 20261002_1550", fontsize="medium")
ax = axes[1]
for c, lab, col in [("gsr_tot", r"Soil-root $\Sigma g_{sr}$", "#8c6d31"), ("gpL2", r"Xylem half-path $2 g_p LAI$", "#08519c"),
                    ("gwL", r"Stem storage $g_{w,stem} LAI$", "#31a354")]:
    ax.semilogy(x, nm[c], "-o", ms=3, color=col, label=lab)
ax.set_ylabel(r"Conductance ($\mu$m s$^{-1}$ MPa$^{-1}$)")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize="small")
ax = axes[2]
for c, lab, col in [("ev", r"Night $E$ (cuticular)", "k"), ("qs_sum", r"Root uptake $\Sigma q_s$", "#8c6d31"),
                    ("qw_stem", r"Stem storage $q_{w,stem}$ (+ = release)", "#31a354")]:
    ax.plot(x, nm[c], "-o", ms=3, color=col, label=lab)
ax.axhline(0, color="0.5", lw=0.6)
ax.set_ylabel(r"Night-mean flux ($\mu$m s$^{-1}$)")
ax.set_xlabel("Night before analysis day")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize="small")
fig.tight_layout()
out = HERE / "night_psi_l_diagnostic.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(out)
