#!/usr/bin/env python
"""
Diagnostic plot of the plant water balance for a lab simulation batch.

Stacked panels share the time axis (burn-in shaded):
    1. Relative soil moisture s per compartment, with shallow and deep VWC / porosity
    2. Water potentials: soil (psi_s), root base (psi_b), xylem (psi_x), leaf (psi_l)
    3. Soil-root conductance gsr per compartment (log scale)
    4. Water fluxes: transpiration, root water uptake and storage fluxes

Storage fluxes are positive when storage discharges into the xylem and negative
when it refills. With flux closure, qs_total = E - qw_stem - qw_leaf.

By default the time series is read from a saved batch. With ``--set``,
``--s-init-scale``, ``--species-set`` or ``--soil-set``, the model is instead re-run
in memory with those parameter overrides, scaled initial soil moisture or species /
soil texture attribute overrides (nothing is saved to results/). Figures go to ``diagnostics/figures/<label>/`` with a ``settings.yaml``
recording the batch or overrides used.

Usage (from the repo root):
    python lab_experiment_simulation/diagnostics/plot_water_balance.py lab_experiment_simulation/config.yaml
    python lab_experiment_simulation/diagnostics/plot_water_balance.py <config> --batch 20261009_0817 --start 2026-08-05 --end 2026-08-08
    python lab_experiment_simulation/diagnostics/plot_water_balance.py <config> --set VWT=1e-6 --label stem_storage_off --sap-flux
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

DIAG_DIR = Path(__file__).resolve().parent
REPO_ROOT = DIAG_DIR.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "sensitivity_analysis"))
sys.path.append(str(DIAG_DIR.parent))  # after sensitivity_analysis/ so its post_processing wins

import soil
import species_traits
from post_processing import _save, set_plot_style  # sensitivity_analysis/post_processing.py
import run_simulation
from run_simulation import simulate

SIDE_COLORS = ["#08519c", "#d95f0e"]


def load_lab_post_processing() -> Any:
    """Import lab_experiment_simulation/post_processing.py (its name clashes with the sensitivity one)."""
    spec = importlib.util.spec_from_file_location("lab_post_processing", DIAG_DIR.parent / "post_processing.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resolve_batch_dir(results_dir: Path, batch: str) -> Path:
    """Return the batch folder for a batch id, or the newest completed batch for ``latest``."""
    if batch == "latest":
        batches = sorted(p for p in results_dir.iterdir() if any((p / "timeseries").glob("*.csv")))
        if not batches:
            raise FileNotFoundError(f"No completed batches in {results_dir}")
        return batches[-1]
    return results_dir / batch


def compartment(df: pd.DataFrame, variable: str, i: int) -> np.ndarray:
    """Values of a per-compartment (JSON list) variable for compartment ``i`` (0-based)."""
    return np.array([json.loads(cell)[i] for cell in df[variable]])


def measured_s(lab: Dict[str, Any], porosity: float) -> pd.DataFrame:
    """Measured relative soil moisture: each side's shallow and deep VWC probes divided by porosity.

    Columns are ``{side}_shallow`` and ``{side}_deep``.
    """
    df = pd.read_csv(REPO_ROOT / lab["vwc_file"], parse_dates=["datetime"]).set_index("datetime").sort_index()
    out = pd.DataFrame(index=df.index)
    for side in lab["sides"]:
        for depth in ("shallow", "deep"):
            out[f"{side}_{depth}"] = df[f"tree{lab['tree']}_{side}_{depth}"] / porosity
    return out


def daily_water_mm(flux: pd.Series) -> pd.Series:
    """Integrate a ground-area flux (µm/s) to mm of water on each calendar day.

    Each sample is taken to cover the interval since the previous sample. A timestamp
    at midnight belongs to the day that just ended.
    """
    flux = flux.dropna().sort_index()
    if len(flux) < 2:
        return pd.Series(dtype=float)
    dt = flux.index.to_series().diff().dt.total_seconds()
    dt.iloc[0] = dt.iloc[1]
    depth_mm = flux.to_numpy() * dt.to_numpy() / 1000.0
    day = (flux.index - pd.Timedelta(seconds=1)).normalize()
    return pd.Series(depth_mm, index=day).groupby(level=0).sum()


def parse_overrides(items: List[str]) -> Dict[str, Any]:
    """``["VWT=1e-6", "root_frac=[0.7, 0.3]"]`` -> ``{"VWT": 1e-6, "root_frac": [0.7, 0.3]}`` (YAML values)."""
    overrides = {}
    for item in items:
        key, _, value = item.partition("=")
        parsed = yaml.safe_load(value)
        if isinstance(parsed, str):
            try:
                parsed = float(parsed)  # PyYAML reads e.g. 1e-6 (no decimal point) as a string
            except ValueError:
                pass
        overrides[key.strip()] = parsed
    return overrides


def scale_initial_s(scale: float) -> None:
    """Multiply the data-derived initial soil moisture used by ``simulate`` by ``scale`` (capped at 1)."""
    data_initial_s = run_simulation.lab_inputs.initial_s

    def scaled(*args: Any, **kwargs: Any) -> List[float]:
        return [min(v * scale, 1.0) for v in data_initial_s(*args, **kwargs)]

    run_simulation.lab_inputs.initial_s = scaled


def plot_water_balance(df: pd.DataFrame, cfg: Dict[str, Any], out_path: Path, subtitle: str = "") -> Path:
    """Draw the four-panel water balance diagnostic for one model time series."""
    sim_cfg, lab = cfg["simulation"], cfg["lab"]
    sides = lab["sides"]
    porosity = getattr(soil, sim_cfg["soil_texture"])().N
    obs = measured_s(lab, porosity)
    obs = obs[(obs.index >= df.index[0]) & (obs.index <= df.index[-1])]
    t = df.index

    fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    ax_s, ax_psi, ax_g, ax_q = axes

    for i, side in enumerate(sides):
        c = SIDE_COLORS[i % len(SIDE_COLORS)]
        ax_s.plot(t, compartment(df, "s", i), color=c, lw=1.2, label=f"Model $s_{i + 1}$ ({side})")
        ax_s.plot(obs.index, obs[f"{side}_shallow"], color=c, lw=0.8, ls=":",
                  label=f"Shallow VWC/$n$ ({side})")
        ax_s.plot(obs.index, obs[f"{side}_deep"], color=c, lw=0.8, ls="--",
                  label=f"Deep VWC/$n$ ({side})")
        ax_psi.plot(t, compartment(df, "psi_s", i), color=c, lw=1.2, label=f"$\\psi_{{s,{i + 1}}}$ ({side})")
        ax_g.plot(t, compartment(df, "gsr", i), color=c, lw=1.2, label=f"$g_{{sr,{i + 1}}}$ ({side})")
        ax_q.plot(t, compartment(df, "qs", i), color=c, lw=1.0, ls="--", label=f"$q_{{s,{i + 1}}}$ ({side})")
    ax_s.set_ylabel("Relative soil\nmoisture $s$ (-)")

    for var, label, color in [("psi_b", "$\\psi_b$ (root base)", "0.4"),
                              ("psi_x", "$\\psi_x$ (xylem)", "#31a354"),
                              ("psi_l", "$\\psi_l$ (leaf)", "#756bb1")]:
        ax_psi.plot(t, df[var], color=color, lw=1.2, label=label)
    ax_psi.set_ylabel("Water potential\n(MPa)")

    ax_g.set_yscale("log")
    ax_g.set_ylabel("Soil-root conductance\n($\\mu$m s$^{-1}$ MPa$^{-1}$)")

    qs_total = np.array([sum(json.loads(cell)) for cell in df["qs"]])
    ax_q.plot(t, qs_total, color="k", lw=1.4, label="$q_s$ total (root uptake)")
    ax_q.plot(t, df["ev"], color="#e31a1c", lw=1.4, label="$E$ (transpiration)")
    ax_q.plot(t, df["qw_stem"], color="#31a354", lw=1.2, label="$q_{w,stem}$ (+ discharge)")
    ax_q.plot(t, df["qw_leaf"], color="#756bb1", lw=1.2, label="$q_{w,leaf}$ (+ discharge)")
    ax_q.axhline(0, color="0.3", lw=0.6)
    ax_q.set_ylabel("Flux, ground area\n($\\mu$m s$^{-1}$)")

    burn_in = df.index[df["is_burn_in"].to_numpy()]
    for ax in axes:
        if len(burn_in):
            ax.axvspan(burn_in[0], burn_in[-1], color="0.9", zorder=0)
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize="small", frameon=False)

    ax_q.xaxis.set_major_locator(mdates.DayLocator())
    ax_q.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax_q.set_xlabel("Date (2026); shaded = burn-in")
    title = f"{sim_cfg['name']}: soil moisture, water potentials and water balance"
    fig.suptitle(f"{title}\n{subtitle}" if subtitle else title)
    return _save(fig, out_path, {"dpi": 200})


def plot_daily_uptake(df: pd.DataFrame, cfg: Dict[str, Any], out_path: Path, subtitle: str = "") -> Path:
    """Bar chart of daily root uptake against daily sap flux (mm of water), one panel per side.

    Burn-in steps are dropped. Both fluxes are ground-area rates integrated over each calendar day.
    Negative root uptake is left out of the daily total; it still appears on the time-series plots.
    """
    lab = cfg["lab"]
    compared = df[~df["is_burn_in"].to_numpy()]
    sap = run_simulation.lab_inputs.sap_flux_ground_basis(lab)
    sap = sap[(sap.index >= compared.index[0]) & (sap.index <= compared.index[-1])]

    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(10, 6), sharex=True)
    for i, (ax, side) in enumerate(zip(np.atleast_1d(axes), lab["sides"])):
        rwu = daily_water_mm(pd.Series(compartment(compared, "qs", i), index=compared.index).clip(lower=0))
        measured = daily_water_mm(sap[side])
        days = rwu.index.union(measured.index).sort_values()
        x = mdates.date2num(days)
        width = 0.35
        color = SIDE_COLORS[i % len(SIDE_COLORS)]
        ax.bar(x - width / 2, rwu.reindex(days), width=width, color=color, label="Model root uptake")
        ax.bar(x + width / 2, measured.reindex(days), width=width, color=color, alpha=0.4, label="Sap flux")
        ax.set_ylabel("Daily total (mm)\n(uptake < 0 excluded)")
        ax.set_title(f"{side} side")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize="small", frameon=False)
    axes[-1].xaxis_date()
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    axes[-1].set_xlabel("Date (2026)")
    title = f"{cfg['simulation']['name']}: daily root uptake vs sap flux"
    fig.suptitle(f"{title}\n{subtitle}" if subtitle else title)
    return _save(fig, out_path, {"dpi": 200})


def main(argv: Optional[List[str]] = None) -> Path:
    """Plot the water balance diagnostic for one batch."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="Path to the lab simulation config.yaml")
    parser.add_argument("--batch", default="latest", help="Batch id under results/, or 'latest'")
    parser.add_argument("--start", help="First date to plot (default: start of run, incl. burn-in)")
    parser.add_argument("--end", help="Last date to plot, inclusive (default: end of run)")
    parser.add_argument("--set", nargs="+", default=[], metavar="PARAM=VALUE",
                        help="Re-run the model in memory with these params overrides instead of reading a batch")
    parser.add_argument("--sap-flux", action="store_true",
                        help="Also save the root uptake vs sap flux figure (post_processing.py style, burn-in dropped)")
    parser.add_argument("--s-init-scale", type=float, default=1.0,
                        help="Re-run in memory with the data-derived initial soil moisture multiplied by this factor")
    parser.add_argument("--species-set", nargs="+", default=[], metavar="ATTR=VALUE",
                        help="Re-run in memory with these species class attributes overridden (e.g. LAI=3)")
    parser.add_argument("--soil-set", nargs="+", default=[], metavar="ATTR=VALUE",
                        help="Re-run in memory with these soil texture class attributes overridden (e.g. KS=10, cm/day)")
    parser.add_argument("--root-area", nargs="+", default=[], metavar="SIDE=CM2",
                        help="Sapwood area per side for the sap flux comparison (e.g. left=2.65 right=1.57)")
    parser.add_argument("--no-hr", action="store_true",
                        help="Clip root uptake at zero so hydraulic redistribution does not rewet a compartment")
    parser.add_argument("--label", help="Figure subfolder under figures/ (default: 'baseline', or the overrides with --set)")
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    set_plot_style(10)

    if args.set or args.s_init_scale != 1.0 or args.species_set or args.soil_set or args.no_hr:
        overrides = parse_overrides(args.set)
        unknown = set(overrides) - set(cfg["params"])
        if unknown:
            raise KeyError(f"Unknown params: {sorted(unknown)}")
        cfg["params"].update(overrides)
        changes = dict(overrides)
        if args.s_init_scale != 1.0:
            scale_initial_s(args.s_init_scale)
            changes["s_init_scale"] = args.s_init_scale
        species_overrides = parse_overrides(args.species_set)
        soil_overrides = parse_overrides(args.soil_set)
        for cls, class_overrides in [(getattr(species_traits, cfg["simulation"]["species"]), species_overrides),
                                     (getattr(soil, cfg["simulation"]["soil_texture"]), soil_overrides)]:
            for attr, value in class_overrides.items():
                if not hasattr(cls, attr):
                    raise KeyError(f"{cls.__name__} has no attribute {attr}")
                setattr(cls, attr, value)
        changes.update(species_overrides)
        changes.update(soil_overrides)
        if args.no_hr:
            cfg["params"]["clip_negative_uptake"] = True
            changes["no_hr"] = True
        print(f"Running in memory with {changes}")
        df = simulate(cfg)
        label = args.label or "_".join(f"{k}={v}" for k, v in changes.items()).replace(" ", "")
        subtitle = "In-memory run: " + ", ".join(f"{k} = {v}" for k, v in changes.items())
        source = {"source": "in-memory run of config", "overrides": overrides, "s_init_scale": args.s_init_scale,
                  "species_overrides": species_overrides, "soil_overrides": soil_overrides,
                  "no_hr": bool(args.no_hr),
                  "root_area_cm2": parse_overrides(args.root_area) or None,
                  "soil_KS_cm_per_day": getattr(soil, cfg["simulation"]["soil_texture"]).KS}
    else:
        batch_dir = resolve_batch_dir(config_path.parent / "results", args.batch)
        df = pd.read_csv(batch_dir / "timeseries" / f"{cfg['simulation']['name']}.csv", parse_dates=["datetime"])
        label = args.label or "baseline"
        subtitle = f"Batch {batch_dir.name}"
        source = {"source": "saved batch", "batch": batch_dir.name}
    if args.root_area:
        cfg["lab"]["root_area_cm2"] = parse_overrides(args.root_area)
    df = df.set_index("datetime")

    out_dir = DIAG_DIR / "figures" / label
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "settings.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump({"config": str(config_path.relative_to(REPO_ROOT)), **source}, f, sort_keys=False)
    if args.start:
        df = df[df.index >= pd.Timestamp(args.start)]
    if args.end:
        df = df[df.index < pd.Timestamp(args.end) + pd.Timedelta(days=1)]

    suffix = f"_{args.start or 'start'}_{args.end or 'end'}" if (args.start or args.end) else ""
    path = plot_water_balance(df, cfg, out_dir / f"water_balance{suffix}.png", subtitle)
    print(f"Saved {path.relative_to(REPO_ROOT)}")
    daily_path = plot_daily_uptake(df, cfg, out_dir / f"daily_uptake_vs_sap_flux{suffix}.png", subtitle)
    print(f"Saved {daily_path.relative_to(REPO_ROOT)}")

    if args.sap_flux:
        lab_pp = load_lab_post_processing()
        plots_cfg = cfg.get("plots", {})
        common = plots_cfg.get("common", {})
        plot_cfg = dict(plots_cfg.get("sap_flux_vs_uptake", {}))
        plot_cfg["title"] = f"{plot_cfg.get('title', '')}\n{subtitle}"
        sap_path = lab_pp.plot_sap_flux_vs_uptake(df[~df["is_burn_in"]], cfg, plot_cfg, common, out_dir)
        if suffix:
            sap_path = sap_path.replace(sap_path.with_name(f"{sap_path.stem}{suffix}{sap_path.suffix}"))
        print(f"Saved {sap_path.relative_to(REPO_ROOT)}")
    return path


if __name__ == "__main__":
    main()
