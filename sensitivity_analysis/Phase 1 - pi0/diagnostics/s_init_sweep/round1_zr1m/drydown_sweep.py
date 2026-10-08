"""s_init sweep, round 1a: S2 med run at s_init 0.30, 0.28, 0.26, 0.25, 0.24, 0.23, both compartments drying.

Ran with 1 m compartments, before the gp lag fix. Showed that s_init <= 0.26 pins psi_l at the
-10 MPa solver bound and that lowering s_init alone gives a pre-stressed start, not a drydown.
"""

import copy
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
import sensitivity_analysis as sa
import post_processing as pp

cfg = sa.load_config(repo / "sensitivity_analysis/Phase 1 - pi0/config.yaml")
sim = dict(cfg["simulation"], zr_arr=[1.0, 1.0])
weather = sa.load_weather(repo / sim["weather_file"], float(sim["timestepM"]),
                          float(sim["burn_in_days"]) + float(sim["analysis_days"]))
run = next(r for r in cfg["runs"] if r["name"] == "med")
s_values = [0.30, 0.28, 0.26, 0.25, 0.24, 0.23]

results = {}
for s0 in s_values:
    params, scenario = sa.resolve_params(cfg, run, "S2")
    params["s_init"] = [s0, s0]
    scenario = copy.deepcopy(scenario)
    scenario["post_burn_soil_dynamics"] = ["drydown", "drydown"]
    model = sa.build_model(params, sim, weather)
    res = sa.run_simulation(model, params, scenario, weather, sim, tag=f"s_init={s0}")
    df = sa.outputs_to_dataframe(res, weather, sim)
    results[s0] = df[~df["is_burn_in"]]

pp.set_plot_style(10)
fig, axes = plt.subplots(4, 1, figsize=(9, 12), sharex=True)
colors = plt.cm.viridis(np.linspace(0, 0.9, len(s_values)))
for c, (s0, df) in zip(colors, results.items()):
    x = df["days_since_burn_in"]
    s1 = np.array([json.loads(v)[0] for v in df["s"]])
    psi1 = np.array([json.loads(v)[0] for v in df["psi_s"]])
    axes[0].plot(x, s1, color=c, label=f"s_init = {s0}")
    axes[1].plot(x, psi1, color=c)
    axes[2].plot(x, df["ev"], color=c)
    axes[3].plot(x, df["psi_l"], color=c)
    daily = df.assign(day=np.floor(x - 1e-9)).groupby("day")
    print(f"s_init={s0}: s end={s1[-1]:.3f}, psi_s end={psi1[-1]:.3f} MPa, "
          f"peak E day1={daily['ev'].max().iloc[0]:.3f} day15={daily['ev'].max().iloc[-1]:.3f}, "
          f"min psi_l day1={daily['psi_l'].min().iloc[0]:.2f} day15={daily['psi_l'].min().iloc[-1]:.2f}")

axes[0].set_ylabel("Relative soil moisture, $s$ (-)")
axes[1].set_ylabel("Soil water potential, $\\psi_s$ (MPa)")
axes[2].set_ylabel("Transpiration, $E$ ($\\mu$m s$^{-1}$)")
axes[3].set_ylabel("Leaf water potential, $\\psi_l$ (MPa)")
axes[3].set_xlabel("Days")
axes[3].set_xlim(0, 15)
axes[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize="small")
fig.suptitle("Drydown sweep: both compartments dry, med pi0 (compartment 1 shown for soil)")
fig.tight_layout()
out = Path(__file__).resolve().parent / "drydown_sweep_s_init.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(out)
