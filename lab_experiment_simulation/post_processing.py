#!/usr/bin/env python
"""
Make figures comparing a lab simulation batch (from run_simulation.py) with lab data.

Plot settings are read from the ``plots`` section of the config at plot time,
so figures can be re-made after editing the config without re-running the model.
To add a figure type, write a ``plot_*`` function with the same signature and
register it in ``PLOT_FUNCTIONS`` under its config key.

Usage (from the repo root):
    python lab_experiment_simulation/post_processing.py lab_experiment_simulation/config.yaml --batch latest
    python lab_experiment_simulation/post_processing.py <config> --batch 20261009_0900 --only sap_flux_vs_uptake
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "sensitivity_analysis"))

from post_processing import _panel_ylim, _save, set_plot_style  # sensitivity_analysis/post_processing.py, first on sys.path

import lab_inputs


def resolve_batch_dir(config_dir: Path, batch: str) -> Path:
    """Return the batch folder for a batch id, or the newest batch for ``latest``."""
    results_dir = config_dir / "results"
    if batch == "latest":
        batches = sorted(p for p in results_dir.iterdir() if any((p / "timeseries").glob("*.csv")))
        if not batches:
            raise FileNotFoundError(f"No completed batches in {results_dir}")
        return batches[-1]
    batch_dir = results_dir / batch
    if not (batch_dir / "timeseries").is_dir():
        raise FileNotFoundError(f"No timeseries folder in {batch_dir}")
    return batch_dir


def load_timeseries(batch_dir: Path, sim_cfg: Dict[str, Any]) -> pd.DataFrame:
    """Model time series for the comparison window (burn-in dropped), indexed by datetime."""
    df = pd.read_csv(batch_dir / "timeseries" / f"{sim_cfg['name']}.csv", parse_dates=["datetime"])
    return df[~df["is_burn_in"]].set_index("datetime")


def _compartment(df: pd.DataFrame, variable: str, i: int) -> np.ndarray:
    """Values of a per-compartment variable for compartment ``i`` (0-based)."""
    return np.array([json.loads(cell)[i] for cell in df[variable]])


def plot_sap_flux_vs_uptake(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                            common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: model root water uptake (qs) per compartment against sap flux on a ground-area
    basis for the matching side, one panel per side. Dotted lines mark irrigation, with salt noted."""
    sim_cfg, lab = cfg["simulation"], cfg["lab"]
    sap = lab_inputs.sap_flux_ground_basis(lab)
    sap = sap[(sap.index >= df.index[0]) & (sap.index <= df.index[-1])]
    irrigation = pd.read_csv(REPO_ROOT / lab["irrigation_file"], parse_dates=["timestamp"]).set_index("timestamp")
    irrigation = irrigation[(irrigation.index >= df.index[0]) & (irrigation.index <= df.index[-1])]
    linewidth = float(common.get("linewidth", 1.2))
    titles = plot_cfg.get("panel_titles") or lab["sides"]

    width, height = common.get("figsize", [9, 3.0])
    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])), sharex=True)
    for i, (ax, side) in enumerate(zip(axes, lab["sides"])):
        ax.plot(df.index, _compartment(df, "qs", i), color=common.get("model_color"), lw=linewidth,
                label=f"Model root uptake $q_{{s,{i + 1}}}$")
        ax.plot(sap.index, sap[side], color=common.get("data_color"), lw=linewidth, marker=".", ms=3,
                label=f"Sap flux, {side} root ($d$ = {lab['root_diameter_cm'][side]} cm)")
        prefix = f"tree{lab['tree']}_{side}"
        for t, row in irrigation[irrigation[f"{prefix}_vol_liters"] > 0].iterrows():
            ax.axvline(t, color="0.4", ls=":", lw=0.8)
            note = f"{row[f'{prefix}_vol_liters']} L"
            if row[f"{prefix}_moles"] > 0:
                note += f", {row[f'{prefix}_moles']} mol"
            ax.annotate(note, (t, 1.0), xycoords=("data", "axes fraction"), fontsize="x-small",
                        ha="left", va="top", rotation=90)
        ax.set_title(titles[i], fontsize="medium")
        ax.set_ylabel(plot_cfg.get("ylabel", "Flux ($\\mu$m s$^{-1}$)"))
        ylim = _panel_ylim(plot_cfg.get("ylim"), i)
        if ylim is not None:
            ax.set_ylim(ylim)
        legend_cfg = common.get("legend", {})
        ax.legend(loc=legend_cfg.get("loc", "upper left"), bbox_to_anchor=legend_cfg.get("bbox_to_anchor", [1.01, 1.0]),
                  fontsize=legend_cfg.get("fontsize", "small"), frameon=legend_cfg.get("frameon", False))

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    axes[-1].set_xlabel("Date (2026)")
    fig.suptitle(plot_cfg.get("title", ""))
    return _save(fig, out_dir / f"sap_flux_vs_uptake.{common.get('format', 'png')}", common)


# Config ``plots`` key -> plotting function.
PLOT_FUNCTIONS: Dict[str, Callable[..., Path]] = {
    "sap_flux_vs_uptake": plot_sap_flux_vs_uptake,
}


def main(argv: Optional[List[str]] = None) -> List[Path]:
    """Make every enabled figure for a batch.

    Returns:
        Paths of the saved figures.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="Path to the lab simulation config.yaml")
    parser.add_argument("--batch", default="latest", help="Batch id under results/, or 'latest'")
    parser.add_argument("--only", nargs="+", help="Subset of plot keys to make")
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    plots_cfg = cfg.get("plots", {})
    common = plots_cfg.get("common", {})
    set_plot_style(common.get("font_size", 10))

    batch_dir = resolve_batch_dir(config_path.parent, args.batch)
    df = load_timeseries(batch_dir, cfg["simulation"])
    print(f"Batch: {batch_dir.name}")

    saved = []
    for key, plot_cfg in plots_cfg.items():
        if key == "common" or not plot_cfg.get("enabled", True):
            continue
        if args.only and key not in args.only:
            continue
        if key not in PLOT_FUNCTIONS:
            raise KeyError(f"No plot function registered for '{key}'")
        path = PLOT_FUNCTIONS[key](df, cfg, plot_cfg, common, batch_dir / "figures")
        print(f"  Saved {path.relative_to(batch_dir)}")
        saved.append(path)
    return saved


if __name__ == "__main__":
    main()
