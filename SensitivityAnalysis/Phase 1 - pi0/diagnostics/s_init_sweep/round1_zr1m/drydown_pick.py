"""s_init sweep, round 1b: S2 at s_init 0.32, 0.30, 0.29 for all pi0 runs, both compartments drying.

Ran with 1 m compartments, before the gp lag fix. Led to s_init = 0.32 (batch 20261002_1521),
later superseded by round 2 with 0.2 m compartments.
"""

import copy
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

repo = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
sys.path.insert(0, str(repo / "SensitivityAnalysis"))
import sensitivity_analysis as sa
import post_processing as pp

cfg = sa.load_config(repo / "SensitivityAnalysis/Phase 1 - pi0/config.yaml")
sim = dict(cfg["simulation"], zr_arr=[1.0, 1.0])
weather = sa.load_weather(repo / sim["weather_file"], float(sim["timestepM"]),
                          float(sim["burn_in_days"]) + float(sim["analysis_days"]))
colors = cfg["plots"]["common"]["run_colors"]
s_values = [0.32, 0.30, 0.29]

pp.set_plot_style(10)
fig, axes = plt.subplots(3, len(s_values), figsize=(5 * len(s_values), 9), sharex=True, sharey="row")
for j, s0 in enumerate(s_values):
    for run in cfg["runs"]:
        params, scenario = sa.resolve_params(cfg, run, "S2")
        params["s_init"] = [s0, s0]
        scenario = copy.deepcopy(scenario)
        scenario["post_burn_soil_dynamics"] = ["drydown", "drydown"]
        model = sa.build_model(params, sim, weather)
        res = sa.run_simulation(model, params, scenario, weather, sim, tag=f"{run['name']} s={s0}")
        df = sa.outputs_to_dataframe(res, weather, sim)
        df = df[~df["is_burn_in"]]
        x = df["days_since_burn_in"]
        psi1 = np.array([json.loads(v)[0] for v in df["psi_s"]])
        axes[0, j].plot(x, psi1, color=colors[run["name"]], label=run["name"])
        axes[1, j].plot(x, df["ev"], color=colors[run["name"]])
        axes[2, j].plot(x, df["psi_l"], color=colors[run["name"]])
        daily = df.assign(day=np.floor(x - 1e-9)).groupby("day")
        print(f"RESULT s_init={s0} {run['name']}: psi_s {psi1[0]:.3f}->{psi1[-1]:.3f}, "
              f"peak E d1={daily['ev'].max().iloc[0]:.3f} d15={daily['ev'].max().iloc[-1]:.3f}, "
              f"min psi_l d1={daily['psi_l'].min().iloc[0]:.2f} d15={daily['psi_l'].min().iloc[-1]:.2f}, "
              f"overall min psi_l={df['psi_l'].min():.2f}")
    axes[0, j].set_title(f"s_init = {s0}")
    axes[2, j].set_xlabel("Days")

axes[0, 0].set_ylabel("Soil water potential,\ncompartment 1 (MPa)")
axes[1, 0].set_ylabel("Transpiration, $E$ ($\\mu$m s$^{-1}$)")
axes[2, 0].set_ylabel("Leaf water potential, $\\psi_l$ (MPa)")
axes[0, -1].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
axes[-1, 0].set_xlim(0, 15)
fig.suptitle("Pre-stressed drydown candidates (both compartments dry)")
fig.tight_layout()
out = Path(__file__).resolve().parent / "drydown_pick_s_init.png"
fig.savefig(out, dpi=130, bbox_inches="tight")
print(out)
