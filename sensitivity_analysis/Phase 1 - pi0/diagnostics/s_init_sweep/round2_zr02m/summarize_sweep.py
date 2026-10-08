"""s_init sweep, round 2: summary table and comparison figure (plus batch 20261002_1550 at 0.32).

Reads the CSVs written by run_sweep.py (0.2 m compartments, before the gp lag fix).
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

REPO = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "sensitivity_analysis"))
from dics import P_ATM, R  # noqa: E402
from post_processing import set_plot_style  # noqa: E402
import species_traits  # noqa: E402

HERE = Path(__file__).resolve().parent
PHASE = REPO / "sensitivity_analysis" / "Phase 1 - pi0"
BATCH = PHASE / "results" / "20261002_1550" / "timeseries"
RUNS = ["low", "med", "high"]
S_INITS = [0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.45]
FIG_S = [float(a) for a in sys.argv[1:]] or [0.32, 0.34, 0.36, 0.38]
GCUT = species_traits.Pecan.GCUT  # mm/s


def load(s_init, run):
    if abs(s_init - 0.32) < 1e-9:
        return pd.read_csv(BATCH / f"{run}__S2.csv")
    return pd.read_csv(HERE / f"s{int(round(s_init * 100))}_{run}.csv")


def first(df, col):
    return np.array([json.loads(c)[0] for c in df[col]])


data = {(s, r): load(s, r) for s in S_INITS for r in RUNS}
records = []
for (s, r), df in data.items():
    an = df[~df["is_burn_in"]].copy()
    an["day"] = np.floor(an["days_since_burn_in"] - 1e-9).astype(int)  # step ending at t belongs to day floor(t-)
    an["gcut_mol"] = GCUT * P_ATM / (1000 * R * an["forcing_ta"])
    daily = an.groupby("day").agg(Epk=("ev", "max"), gsw_pk=("gsw", "max"), gcut=("gcut_mol", "max"),
                                  psil_min=("psi_l", "min"))
    floor_days = daily.index[daily["gsw_pk"] < 1.05 * daily["gcut"]]
    psi_s1 = first(an, "psi_s")
    records.append({
        "s_init": s, "run": r,
        "s_final": first(an, "s")[-1], "psi_s_final": psi_s1[-1],
        "Epk_d1": daily["Epk"].iloc[0], "Epk_d15": daily["Epk"].iloc[-1],
        "floor_day_gsw": int(floor_days[0]) + 1 if len(floor_days) else np.nan,
        "floor_day_E": (int(np.flatnonzero(daily["Epk"].to_numpy() < 1.03 * 0.0251)[0]) + 1
                        if (daily["Epk"] < 1.03 * 0.0251).any() else np.nan),
        "psil_min": an["psi_l"].min(),
        "psil_predawn_d15": an[np.isclose((an["time_days"] % 1) * 24, 5.0)]["psi_l"].iloc[-1],
        "n_psil_floor": int((df["psi_l"] <= -9.99).sum()),
        "max_abs_fbal": df["flux_balance"].abs().max(),
        "nan": int(df[["psi_l", "ev"]].isna().sum().sum()),
        "_daily_Epk": daily["Epk"].to_numpy(), "_daily_psil": daily["psil_min"].to_numpy(),
    })
tab = pd.DataFrame(records)

# separation between pi0 runs: daily peak E spread (high-low)/med, and daily min psi_l spread
sep = []
for s in S_INITS:
    sub = tab[tab["s_init"] == s].set_index("run")
    e = np.vstack(sub.loc[RUNS, "_daily_Epk"])
    p = np.vstack(sub.loc[RUNS, "_daily_psil"])
    e_spread = (e.max(0) - e.min(0)) / e[1]
    merged = np.flatnonzero(e_spread < 0.05)
    sep.append({"s_init": s,
                "E_spread_d1_%": 100 * e_spread[0], "E_spread_d8_%": 100 * e_spread[7],
                "E_spread_d15_%": 100 * e_spread[-1],
                "E_merge_day": int(merged[0]) + 1 if len(merged) else np.nan,
                "psil_spread_d1": p[0, 0] - p[2, 0], "psil_spread_d15": p[0, -1] - p[2, -1],
                "psil_low_min": p[0].min(), "psil_med_min": p[1].min(), "psil_high_min": p[2].min()})
sep = pd.DataFrame(sep)
pd.set_option("display.width", 300, "display.precision", 4, "display.max_columns", 30)
print(tab.drop(columns=["_daily_Epk", "_daily_psil"]).to_string(index=False))
print()
print(sep.to_string(index=False))
print("\nDaily peak E (um/s) per s_init/run:")
for s in S_INITS:
    for r in RUNS:
        v = tab[(tab.s_init == s) & (tab.run == r)]["_daily_Epk"].iloc[0]
        print(f"  {s} {r:4s}", " ".join(f"{x:.3f}" for x in v))

# figure
with open(PHASE / "config.yaml", encoding="utf-8") as f:
    common = yaml.safe_load(f)["plots"]["common"]
set_plot_style(common.get("font_size", 10))
colors = common["run_colors"]
fig, axes = plt.subplots(3, len(FIG_S), figsize=(4.2 * len(FIG_S), 8.5), sharex=True, sharey="row")
for j, s in enumerate(FIG_S):
    for r in RUNS:
        df = data[(s, r)]
        an = df[~df["is_burn_in"]]
        x = an["days_since_burn_in"]
        axes[0, j].plot(x, first(an, "psi_s"), color=colors[r], lw=1.2, label=r)
        axes[1, j].plot(x, an["ev"], color=colors[r], lw=0.9, label=r)
        axes[2, j].plot(x, an["psi_l"], color=colors[r], lw=0.9, label=r)
    tag = " (batch 20261002_1550)" if abs(s - 0.32) < 1e-9 else ""
    axes[0, j].set_title(f"$s_{{init}}$ = {s:.2f}{tag}", fontsize="medium")
    axes[2, j].set_xlabel("Days since burn-in")
    axes[2, j].set_xlim(0, 15)
axes[0, 0].set_ylabel("Soil $\\psi_s$, comp. 1 (MPa)")
axes[1, 0].set_ylabel("Transpiration $E$ ($\\mu$m s$^{-1}$)")
axes[2, 0].set_ylabel("Leaf $\\psi_l$ (MPa)")
handles, _ = axes[1, -1].get_legend_handles_labels()
labs = ["low: $\\pi_{0,leaf}$ -1.1, $\\pi_{0,stem}$ -0.9", "med: -1.4 / -1.1", "high: -1.7 / -1.3"]
axes[0, -1].legend(handles, labs, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize="small", frameon=False)
fig.suptitle("S2 drydown, zr = species (0.2 m): initial soil moisture sweep")
fig.tight_layout()
out = HERE / "sinit_sweep_comparison.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(out)
