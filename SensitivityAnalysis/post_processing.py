#!/usr/bin/env python
"""
Make figures from a sensitivity analysis batch produced by sensitivity_analysis.py.

Plot settings (axis limits, labels, legend, colors, window) are read from the
``plots`` section of the Phase config at plot time, so figures can be
re-made after editing the config without re-running the model.

Each figure type has its own ``plot_*`` function. To add a figure type,
write a new function with the same signature and register it in
``PLOT_FUNCTIONS`` under the key used in the config's ``plots`` section.

Usage (from the repo root):
    python SensitivityAnalysis/post_processing.py "SensitivityAnalysis/Phase 1 - pi0/config.yaml" --batch latest
    python SensitivityAnalysis/post_processing.py <config> --batch 20261001_1200 --only psi_l
"""

import argparse
import json
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

SeriesDict = Dict[Tuple[str, str], pd.DataFrame]


def set_plot_style(font_size: float = 10) -> None:
    """Apply the project plot style (seaborn whitegrid, serif font stack)."""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Caladea", "STIX Two Text", "Constantia", "Georgia", "Cambria"],
        "font.size": font_size,
        "mathtext.fontset": "stix",
    })


def run_label(overrides_json: str, display_names: Dict[str, str]) -> str:
    """Legend label from a run's overrides, e.g. ``$\\pi_{0,leaf}$ = -1.1; $\\pi_{0,stem}$ = -0.9``.

    Args:
        overrides_json: The manifest ``overrides`` cell (JSON dict of parameter: value).
        display_names: Parameter name -> display string (``plots.common.display_names``).
    """
    overrides = json.loads(overrides_json)
    return "; ".join(f"{display_names.get(k, k)} = {v}" for k, v in overrides.items())


def resolve_batch_dir(phase_dir: Path, batch: str) -> Path:
    """Return the batch folder for a batch id, or the newest batch for ``latest``."""
    results_dir = phase_dir / "results"
    if batch == "latest":
        batches = sorted(p for p in results_dir.iterdir() if (p / "manifest.csv").exists())
        if not batches:
            raise FileNotFoundError(f"No completed batches in {results_dir}")
        return batches[-1]
    batch_dir = results_dir / batch
    if not (batch_dir / "manifest.csv").exists():
        raise FileNotFoundError(f"No manifest.csv in {batch_dir}")
    return batch_dir


def load_batch(batch_dir: Path) -> Tuple[pd.DataFrame, SeriesDict]:
    """Load the manifest and every time series in a batch.

    Args:
        batch_dir: A ``results/<batch_id>/`` folder.

    Returns:
        ``(manifest, series)`` where ``series`` maps ``(run, scenario)`` to a DataFrame.
    """
    manifest = pd.read_csv(batch_dir / "manifest.csv")
    series = {
        (row.run, row.scenario): pd.read_csv(batch_dir / row.timeseries_file)
        for row in manifest.itertuples()
    }
    return manifest, series


def _window(df: pd.DataFrame, window_days: Optional[List[float]]) -> pd.DataFrame:
    """Rows of ``df`` inside ``window_days`` (days since end of burn-in)."""
    if window_days is None:
        return df[~df["is_burn_in"]]
    lo, hi = window_days
    x = df["days_since_burn_in"]
    return df[(x >= lo) & (x <= hi)]


def _values(df: pd.DataFrame, variable: str, compartment: Optional[int]) -> np.ndarray:
    """Values of ``variable``; for per-compartment pairs, pick ``compartment`` (1-based)."""
    if variable not in df.columns:
        raise KeyError(f"Variable '{variable}' not in time series")
    if compartment is None:
        return df[variable].to_numpy(dtype=float)
    return np.array([json.loads(cell)[compartment - 1] for cell in df[variable]])


def _panel_ylim(ylim: Any, panel_index: int) -> Optional[List[float]]:
    """Y limits for one panel; ``ylim`` is null, one ``[min, max]`` pair, or a list of pairs per panel."""
    if ylim is None:
        return None
    if any(isinstance(item, (list, tuple)) for item in ylim):
        return ylim[panel_index]
    return ylim


def _plot_timeseries(
    manifest: pd.DataFrame,
    series: SeriesDict,
    scenario: str,
    plot_cfg: Dict[str, Any],
    common: Dict[str, Any],
    out_path: Path,
    panels: List[Dict[str, Any]],
) -> Path:
    """Shared helper: stacked panels with one line per run, saved to ``out_path``.

    Args:
        manifest: Batch manifest (``overrides`` gives the legend entries).
        series: Time series keyed by ``(run, scenario)``.
        scenario: Scenario to plot.
        plot_cfg: The figure's block from the config ``plots`` section.
        common: The config ``plots.common`` block.
        out_path: File to write.
        panels: One dict per panel with ``variable``, optional ``compartment``
            (1-based, for per-compartment pairs) and optional ``title``.

    Returns:
        Path of the saved figure.
    """
    rows = manifest[manifest["scenario"] == scenario]
    display_names = common.get("display_names", {})
    colors = common.get("run_colors", {})
    scale = float(plot_cfg.get("scale", 1.0))
    linewidth = float(common.get("linewidth", 1.2))
    legend_cfg = {**common.get("legend", {}), **plot_cfg.get("legend", {})}

    width, height = common.get("figsize", [8, 3.5])
    fig, axes = plt.subplots(len(panels), 1, figsize=(width, height * len(panels)), sharex=True, squeeze=False)
    axes = axes[:, 0]

    for p, (ax, panel) in enumerate(zip(axes, panels)):
        for row in rows.itertuples():
            df = _window(series[(row.run, scenario)], common.get("window_days"))
            ax.plot(
                df["days_since_burn_in"],
                _values(df, panel["variable"], panel.get("compartment")) * scale,
                color=colors.get(row.run), linewidth=linewidth,
                label=run_label(row.overrides, display_names),
            )
        ax.set_ylabel(plot_cfg.get("ylabel", panel["variable"]))
        ylim = _panel_ylim(plot_cfg.get("ylim"), p)
        if ylim is not None:
            ax.set_ylim(ylim)
        if panel.get("title"):
            ax.set_title(panel["title"], fontsize="medium")

    axes[0].legend(
        title=legend_cfg.get("title"),
        loc=legend_cfg.get("loc", "upper left"),
        bbox_to_anchor=legend_cfg.get("bbox_to_anchor", [1.01, 1.0]),
        fontsize=legend_cfg.get("fontsize", "small"),
        frameon=legend_cfg.get("frameon", False),
    )

    if common.get("window_days") is not None:
        axes[-1].set_xlim(common["window_days"])
    axes[-1].set_xlabel(common.get("xlabel", "Days since end of burn-in"))
    fig.suptitle(f"{plot_cfg.get('title', '')} ({rows['scenario_label'].iloc[0]})")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=common.get("dpi", 300), bbox_inches="tight")
    plt.close(fig)
    return out_path


def _panels_from_config(plot_cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build panel specs from ``variables``/``compartments``/``panel_titles`` in a figure block."""
    variables = plot_cfg["variables"]
    compartments = plot_cfg.get("compartments")
    if compartments:
        specs = [{"variable": variables[0], "compartment": c} for c in compartments]
    else:
        specs = [{"variable": v} for v in variables]
    for spec, title in zip(specs, plot_cfg.get("panel_titles") or []):
        spec["title"] = title
    return specs


def plot_transpiration(manifest: pd.DataFrame, series: SeriesDict, scenario: str,
                       plot_cfg: Dict[str, Any], common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: time series of transpiration (ev), one line per run, one figure per scenario."""
    return _plot_timeseries(
        manifest, series, scenario, plot_cfg, common,
        out_dir / f"transpiration.{common.get('format', 'png')}",
        panels=_panels_from_config(plot_cfg),
    )


def plot_root_uptake(manifest: pd.DataFrame, series: SeriesDict, scenario: str,
                     plot_cfg: Dict[str, Any], common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: stacked panels of root water uptake (qs), one panel per soil compartment, one line per run; one figure per scenario."""
    return _plot_timeseries(
        manifest, series, scenario, plot_cfg, common,
        out_dir / f"root_uptake.{common.get('format', 'png')}",
        panels=_panels_from_config(plot_cfg),
    )


def plot_soil_water_potential(manifest: pd.DataFrame, series: SeriesDict, scenario: str,
                              plot_cfg: Dict[str, Any], common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: stacked panels of soil water potential (psi_s), one panel per soil compartment, one line per run; one figure per scenario."""
    return _plot_timeseries(
        manifest, series, scenario, plot_cfg, common,
        out_dir / f"soil_water_potential.{common.get('format', 'png')}",
        panels=_panels_from_config(plot_cfg),
    )


def plot_storage_fluxes(manifest: pd.DataFrame, series: SeriesDict, scenario: str,
                        plot_cfg: Dict[str, Any], common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: stacked panels of storage-to-xylem fluxes (qw_stem, qw_leaf), one panel per variable, one line per run; one figure per scenario."""
    return _plot_timeseries(
        manifest, series, scenario, plot_cfg, common,
        out_dir / f"storage_fluxes.{common.get('format', 'png')}",
        panels=_panels_from_config(plot_cfg),
    )


def plot_psi_l(manifest: pd.DataFrame, series: SeriesDict, scenario: str,
               plot_cfg: Dict[str, Any], common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: time series of leaf water potential (psi_l), one line per run, one figure per scenario."""
    return _plot_timeseries(
        manifest, series, scenario, plot_cfg, common,
        out_dir / f"psi_l.{common.get('format', 'png')}",
        panels=_panels_from_config(plot_cfg),
    )


# Config ``plots`` key -> plotting function.
PLOT_FUNCTIONS: Dict[str, Callable[..., Path]] = {
    "transpiration": plot_transpiration,
    "root_uptake": plot_root_uptake,
    "soil_water_potential": plot_soil_water_potential,
    "storage_fluxes": plot_storage_fluxes,
    "psi_l": plot_psi_l,
}


def main(argv: Optional[List[str]] = None) -> List[Path]:
    """Make every enabled figure for every scenario in a batch.

    Returns:
        Paths of the saved figures.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="Path to the Phase config.yaml")
    parser.add_argument("--batch", default="latest", help="Batch id under results/, or 'latest'")
    parser.add_argument("--only", nargs="+", help="Subset of plot keys to make")
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    with open(config_path, "r", encoding="utf-8") as f:
        plots_cfg = yaml.safe_load(f).get("plots", {})
    common = plots_cfg.get("common", {})
    set_plot_style(common.get("font_size", 10))

    batch_dir = resolve_batch_dir(config_path.parent, args.batch)
    manifest, series = load_batch(batch_dir)
    print(f"Batch: {batch_dir.name} ({len(series)} time series)")

    saved = []
    for key, plot_cfg in plots_cfg.items():
        if key == "common" or not plot_cfg.get("enabled", True):
            continue
        if args.only and key not in args.only:
            continue
        if key not in PLOT_FUNCTIONS:
            raise KeyError(f"No plot function registered for '{key}'")
        for scenario in manifest["scenario"].unique():
            path = PLOT_FUNCTIONS[key](manifest, series, scenario, plot_cfg, common, batch_dir / "figures" / scenario)
            print(f"  Saved {path.relative_to(batch_dir)}")
            saved.append(path)
    return saved


if __name__ == "__main__":
    main()
