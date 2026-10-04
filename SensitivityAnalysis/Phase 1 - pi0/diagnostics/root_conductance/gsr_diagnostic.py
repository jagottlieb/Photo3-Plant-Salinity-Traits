"""Root conductance diagnostics.

Part A: soil-root conductance (gsr) vs soil moisture and soil water potential, using the
species root depth (Pecan.ZR = 0.2 m) as hydraulics.py does; compared with plant gp*LAI.
Part B: steps through S2 at s_init = 0.23 with the original 1 m compartments, before the gp
lag fix, to show the lead-up to psi_l pinning at the -10 MPa solver bound. Rerunning Part B
now uses the fixed solver, so it will not exactly reproduce psi_l_collapse_diagnostic.png.
"""

import copy
import json
import sys
from math import pi
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

repo = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(repo / "SensitivityAnalysis"))
import sensitivity_analysis as sa
import post_processing as pp
from soil import ConstantSoil, DrydownSoil

cfg = sa.load_config(repo / "SensitivityAnalysis/Phase 1 - pi0/config.yaml")
sim = cfg["simulation"]
weather = sa.load_weather(repo / sim["weather_file"], float(sim["timestepM"]),
                          float(sim["burn_in_days"]) + float(sim["analysis_days"]))
run = next(r for r in cfg["runs"] if r["name"] == "med")
params, scenario = sa.resolve_params(cfg, run, "S2")
model = sa.build_model(params, sim, weather)
hydro, soil, sp = model["hydro"], model["plant"].soil, model["species"]
zr = hydro.zr
B, rf = params["B_param"], np.array(params["root_frac"], dtype=float)
lai, gpmax = hydro.lai, sp.GPMAX
print(f"LAI={lai}, GPMAX={gpmax:.4f} um/(MPa s) per leaf, plant gp*LAI={gpmax*lai:.3f}, "
      f"GCUT={sp.GCUT:.4f} mm/s, PSILA0={model['plant'].photo.PSILA0:.3f}, PSILA1={model['plant'].photo.PSILA1}")

# ---------------- Part A: gsr vs soil moisture ----------------
pp.set_plot_style(10)
s_grid = np.linspace(0.195, 0.6, 400)
gsr_tot = np.array([hydro.gsr(soil, [s, s], zr, B, rf).sum() for s in s_grid])
rr, kr = 0.2e-3, 1e-8
kr_limit = sum(kr * 101.9e6 * 2 * pi * rr * B * f * zr for f in rf)
psi_s = np.array([soil.stype[0].psi_s(s) for s in s_grid])

fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, x, xl in [(axes[0], s_grid, "Relative soil moisture, $s$ (-)"),
                  (axes[1], -psi_s, "Soil water potential, $-\\psi_s$ (MPa)")]:
    ax.semilogy(x, gsr_tot, color="#08519c", lw=1.8, label="Soil-root, $g_{sr}$ (both compartments)")
    ax.axhline(kr_limit, color="0.5", ls=":", label="Root radial limit (wet soil)")
    for psil, ls in [(0.0, "--"), (-1.5, "-."), (-3.0, (0, (1, 3)))]:
        ax.axhline(gpmax * lai * np.exp(-(psil / 2) ** 2), color="#d94801", ls=ls,
                   label=f"Plant xylem, $g_p\\cdot$LAI at $\\psi_l$ = {psil}")
    ax.set_xlabel(xl)
for s0, lab in [(0.45, "0.45"), (0.32, "0.32 (final S2)"), (0.23, "0.23 (fails)")]:
    axes[0].axvline(s0, color="0.3", lw=0.8)
    axes[0].text(s0, axes[0].get_ylim()[0] if False else 1e-6, f" {lab}", rotation=90, va="bottom", fontsize=8)
    axes[1].axvline(-soil.stype[0].psi_s(s0), color="0.3", lw=0.8)
axes[1].set_xscale("log")
axes[0].set_ylabel("Conductance ($\\mu$m s$^{-1}$ MPa$^{-1}$, ground area)")
axes[1].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize="small")
fig.suptitle(f"Soil-root conductance vs soil water (Berger soil, root depth ZR = {zr} m, B = 10000)")
fig.tight_layout()
fig.savefig(HERE / "gsr_vs_soil_water.png", dpi=140, bbox_inches="tight")
for s0 in [0.45, 0.32, 0.30, 0.26, 0.23, 0.21]:
    print(f"s={s0}: gsr_tot={hydro.gsr(soil, [s0, s0], zr, B, rf).sum():.4g}, psi_s={soil.stype[0].psi_s(s0):.3f}")

for s0 in np.arange(0.36, 0.27, -0.01):
    print(f"s={s0:.2f}: gsr_tot={hydro.gsr(soil, [s0, s0], zr, B, rf).sum():.4g}")

# ---------------- Part B: step through the failing run ----------------
sim = dict(sim, zr_arr=[1.0, 1.0])
params["s_init"] = [0.23, 0.23]
model = sa.build_model(params, sim, weather)
plant, hydro = model["plant"], model["hydro"]
dt = sim["timestepM"] * 60.0
burn = int(sa.steps(float(sim["burn_in_days"]), int(sim["timestepM"])))
n = int(sa.steps(float(sim["burn_in_days"]) + float(sim["analysis_days"]), int(sim["timestepM"])))
plant.soil.dynamics = [ConstantSoil(), ConstantSoil()]
cost, gsr_t = [], []
for i in range(n):
    if i == burn:
        plant.soil.dynamics = [DrydownSoil(), DrydownSoil()]
    plant.update(dt, weather["phi"][i], weather["ta"][i], weather["qa"][i])
    cost.append(hydro.out.cost)
out = plant.output()
t = (np.arange(n) + 1) * dt / 86400 - burn * dt / 86400
psi_l = np.array(out["psi_l"]); ev = np.array(out["ev"]); gp = np.array(out["gp"])
gsw = np.array(out["gsw"]); gsr = np.array(out["gsr"]).sum(axis=0)
qs_arr = np.array(out["qs"]); qs = qs_arr.sum(axis=1) if qs_arr.shape[0] == n else qs_arr.sum(axis=0); qw = np.array(out["qw_stem"]); wst = np.array(out["w_stem"])
psi_b = np.array(out["psi_b"]); psi_x = np.array(out["psi_x"]); cost = np.array(cost)
gcut_mol = sp.GCUT * 101325 / (1000 * 8.314 * np.array(weather["ta"][:n]))

k = np.argmax((psi_l < -9.9) & (t > 0))
print(f"\nCollapse at step {k}, day {t[k]:.2f} post burn-in")
for j in range(k - 6, k + 3):
    print(f"d={t[j]:6.2f} psi_l={psi_l[j]:7.2f} psi_x={psi_x[j]:7.2f} psi_b={psi_b[j]:7.2f} "
          f"E={ev[j]:.4f} qs={qs[j]:+.4f} qw_stem={qw[j]:+.4f} w_stem={wst[j]:.3f} "
          f"gsw={gsw[j]:.4f} (cut={gcut_mol[j]:.4f}) gsr={gsr[j]:.4g} gp*LAI={gp[j]*lai:.4g} cost={cost[j]:.2e}")

w = (t > t[k] - 3) & (t < t[k] + 1)
fig, ax = plt.subplots(5, 1, figsize=(9, 12), sharex=True)
ax[0].plot(t[w], psi_l[w], label="$\\psi_l$"); ax[0].plot(t[w], psi_x[w], label="$\\psi_x$")
ax[0].plot(t[w], psi_b[w], label="$\\psi_b$ (root base)")
ax[0].axhline(model["plant"].photo.PSILA0, color="0.4", ls=":", label="PSILA0 (stomata fully shut)")
ax[0].set_ylabel("MPa"); ax[0].set_ylim(-11, 0)
ax[1].plot(t[w], gsw[w], label="$g_{sw}$ total"); ax[1].plot(t[w], gcut_mol[w], ls="--", label="cuticular part")
ax[1].set_ylabel("mol m$^{-2}$ s$^{-1}$")
ax[2].semilogy(t[w], gsr[w], label="$g_{sr}$ total"); ax[2].semilogy(t[w], gp[w] * lai, label="$g_p\\cdot$LAI")
ax[2].set_ylabel("Conductance")
ax[3].plot(t[w], ev[w], label="E (demand)"); ax[3].plot(t[w], qs[w], label="$q_s$ (root uptake)")
ax[3].plot(t[w], qw[w], label="$q_{w,stem}$ (storage)"); ax[3].set_ylabel("$\\mu$m s$^{-1}$")
ax[4].plot(t[w], wst[w], label="stem storage fraction"); ax[4].semilogy() if False else None
ax4b = ax[4].twinx(); ax4b.semilogy(t[w], np.maximum(cost[w], 1e-30), color="C3", lw=0.8)
ax4b.set_ylabel("solver residual (cost)", color="C3")
ax[4].set_xlabel("Days")
for a in ax:
    a.legend(loc="upper left", bbox_to_anchor=(1.08, 1.0), frameon=False, fontsize="small")
fig.suptitle("s_init = 0.23, med pi0: lead-up to psi_l collapse")
fig.tight_layout()
fig.savefig(HERE / "psi_l_collapse_diagnostic.png", dpi=130, bbox_inches="tight")
print("saved")
