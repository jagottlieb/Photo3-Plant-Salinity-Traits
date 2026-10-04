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
    axes[-1].set_xlabel(common.get("xlabel", "Days"))
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


def _shade(ax: plt.Axes, x: np.ndarray, mask: np.ndarray) -> None:
    """Grey spans where ``mask`` is True (contiguous runs of timesteps)."""
    dt = np.median(np.diff(x))
    edges = np.flatnonzero(np.diff(np.r_[0, mask.astype(int), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        ax.axvspan(x[a] - dt / 2, x[b - 1] + dt / 2, color="0.9", lw=0, zorder=0)


def plot_flux_partition(manifest: pd.DataFrame, series: SeriesDict, scenario: str,
                        plot_cfg: Dict[str, Any], common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: transpiration supply partition, one column per run; one figure per scenario.

    Rows: (1) absolute fluxes: E, root uptake (sum of qs over compartments), qw_stem
    and qw_leaf; (2) per-timestep signed shares of E (positive = supplies E,
    negative = storage refilling, so the root share can exceed 100%), nights
    shaded; (3) daily totals of each flux as % of daily E. A thin black line in
    rows 2-3 is the sum of the shares, a check that the water balance closes (100%).

    Optional ``plot_cfg`` keys: ``ylabels`` (one per row), ``flux_ylim``,
    ``share_ylim``, ``daily_ylim``, ``component_colors`` and ``component_labels``
    (keys root, stem, leaf), ``show_sum``.
    """
    rows = manifest[manifest["scenario"] == scenario]
    display_names = common.get("display_names", {})
    run_colors = common.get("run_colors", {})
    linewidth = float(common.get("linewidth", 1.2))
    legend_cfg = {**common.get("legend", {}), **plot_cfg.get("legend", {})}
    comp_colors = {"root": "#8c6d31", "stem": "#31a354", "leaf": "#a1d99b", **plot_cfg.get("component_colors", {})}
    comp_labels = {"root": "Root uptake $\\Sigma q_s$", "stem": "Stem storage $q_{w,stem}$",
                   "leaf": "Leaf storage $q_{w,leaf}$", **plot_cfg.get("component_labels", {})}
    ylabels = plot_cfg.get("ylabels") or ["Flux ($\\mu$m s$^{-1}$)", "Share of $E$ per step (%)", "Share of daily $E$ (%)"]
    show_sum = plot_cfg.get("show_sum", True)

    width, height = common.get("figsize", [8, 3.5])
    fig, axes = plt.subplots(3, len(rows), figsize=(0.75 * width * len(rows), 0.8 * height * 3),
                             sharex="col", squeeze=False)

    for j, row in enumerate(rows.itertuples()):
        df = _window(series[(row.run, scenario)], common.get("window_days"))
        x = df["days_since_burn_in"].to_numpy()
        ev = df["ev"].to_numpy(dtype=float)
        n_comp = len(json.loads(df["qs"].iloc[0]))
        flux = {
            "root": sum(_values(df, "qs", c) for c in range(1, n_comp + 1)),
            "stem": df["qw_stem"].to_numpy(dtype=float),
            "leaf": df["qw_leaf"].to_numpy(dtype=float),
        }

        ax = axes[0, j]
        ax.plot(x, ev, color="k", lw=linewidth, label="Transpiration $E$")
        for k, q in flux.items():
            ax.plot(x, q, color=comp_colors[k], lw=linewidth, label=comp_labels[k])
        ax.axhline(0, color="0.5", lw=0.6)
        if plot_cfg.get("flux_ylim"):
            ax.set_ylim(plot_cfg["flux_ylim"])
        ax.set_title(run_label(row.overrides, display_names), fontsize="medium", color=run_colors.get(row.run, "k"))

        ax = axes[1, j]
        share = {k: 100 * q / ev for k, q in flux.items()}
        pos_base, neg_base = np.zeros_like(ev), np.zeros_like(ev)
        for k, sh in share.items():
            pos, neg = np.clip(sh, 0, None), np.clip(sh, None, 0)
            ax.fill_between(x, pos_base, pos_base + pos, color=comp_colors[k], lw=0, step="mid", label=comp_labels[k])
            ax.fill_between(x, neg_base, neg_base + neg, color=comp_colors[k], lw=0, step="mid")
            pos_base, neg_base = pos_base + pos, neg_base + neg
        if show_sum:
            ax.plot(x, sum(share.values()), color="k", lw=0.6, label="Sum of shares")
        ax.axhline(100, color="k", lw=0.6, ls=":")
        ax.axhline(0, color="0.5", lw=0.6)
        _shade(ax, x, df["forcing_phi"].to_numpy(dtype=float) <= 0)
        if plot_cfg.get("share_ylim"):
            ax.set_ylim(plot_cfg["share_ylim"])

        ax = axes[2, j]
        dt = np.median(np.diff(x))
        day_idx = np.floor(x - dt / 2 + 1e-9).astype(int)
        days, counts = np.unique(day_idx, return_counts=True)
        days = days[counts == counts.max()]
        e_day = np.array([ev[day_idx == d].sum() for d in days])
        pos_base, neg_base, total = np.zeros(len(days)), np.zeros(len(days)), np.zeros(len(days))
        for k, q in flux.items():
            sh = 100 * np.array([q[day_idx == d].sum() for d in days]) / e_day
            pos, neg = np.clip(sh, 0, None), np.clip(sh, None, 0)
            ax.bar(days + 0.5, pos, bottom=pos_base, width=0.8, color=comp_colors[k])
            ax.bar(days + 0.5, neg, bottom=neg_base, width=0.8, color=comp_colors[k])
            pos_base, neg_base, total = pos_base + pos, neg_base + neg, total + sh
        if show_sum:
            ax.plot(days + 0.5, total, "k.-", lw=0.6, ms=3)
        ax.set_ylim(plot_cfg.get("daily_ylim") or [min(neg_base.min(), 0) - 10, max(pos_base.max(), total.max()) + 10])
        ax.axhline(100, color="k", lw=0.6, ls=":")
        ax.axhline(0, color="0.5", lw=0.6)

    for r, label in enumerate(ylabels):
        axes[r, 0].set_ylabel(label)
    for ax in axes[-1]:
        ax.set_xlabel(common.get("xlabel", "Days"))
        if common.get("window_days") is not None:
            ax.set_xlim(common["window_days"])
    legend_kw = dict(loc=legend_cfg.get("loc", "upper left"), bbox_to_anchor=legend_cfg.get("bbox_to_anchor", [1.01, 1.0]),
                     fontsize=legend_cfg.get("fontsize", "small"), frameon=legend_cfg.get("frameon", False))
    axes[0, -1].legend(**legend_kw)
    handles, labels = axes[1, -1].get_legend_handles_labels()
    axes[1, -1].legend(handles + [plt.Rectangle((0, 0), 1, 1, color="0.9")], labels + ["Night"], **legend_kw)
    fig.suptitle(f"{plot_cfg.get('title', '')} ({rows['scenario_label'].iloc[0]})")
    fig.tight_layout()

    out_path = out_dir / f"flux_partition.{common.get('format', 'png')}"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=common.get("dpi", 300), bbox_inches="tight")
    plt.close(fig)
    return out_path


# Config ``plots`` key -> plotting function.
PLOT_FUNCTIONS: Dict[str, Callable[..., Path]] = {
    "transpiration": plot_transpiration,
    "root_uptake": plot_root_uptake,
    "soil_water_potential": plot_soil_water_potential,
    "storage_fluxes": plot_storage_fluxes,
    "psi_l": plot_psi_l,
    "flux_partition": plot_flux_partition,
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
