#!/usr/bin/env python
"""
Make figures comparing a lab simulation batch (from run_simulation.py) with lab data.

Plot settings are read from the ``plots`` section of the config at plot time,
so figures can be re-made after editing the config without re-running the model.
To add a figure type, write a ``plot_*`` function with the same signature and
register it in ``PLOT_FUNCTIONS`` under its config key.

Usage (from the repo root):
    python lab_experiment_simulation/post_processing.py lab_experiment_simulation/drydown/config.yaml --batch latest
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
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Rectangle
from matplotlib.text import Text as MplText
import numpy as np
import pandas as pd
import yaml

# Deep TEROS probes begin failing ~16:00 on 8/17 and are gone by 17:00 (both sides).
DEEP_PROBE_CRASH = pd.Timestamp("2026-08-17 16:00:00")

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "sensitivity_analysis"))
sys.path.append(str(REPO_ROOT))

from post_processing import _panel_ylim, _save, set_plot_style  # sensitivity_analysis/post_processing.py, first on sys.path

import lab_inputs
import soil


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


def _legend(ax: plt.Axes, common: Dict[str, Any],
            extra_handles: Optional[List[Any]] = None,
            extra_labels: Optional[List[str]] = None,
            handler_map: Optional[Dict[Any, Any]] = None) -> None:
    """Place the panel legend from ``plots.common.legend``."""
    legend_cfg = common.get("legend", {})
    handles, labels = ax.get_legend_handles_labels()
    kwargs: Dict[str, Any] = {
        "loc": legend_cfg.get("loc", "upper left"),
        "bbox_to_anchor": legend_cfg.get("bbox_to_anchor", [1.01, 1.0]),
        "fontsize": legend_cfg.get("fontsize", "small"),
        "frameon": legend_cfg.get("frameon", False),
        "handler_map": handler_map or {},
    }
    if extra_handles:
        # Blank spacer row, then the pulse box dropped slightly below the series entries.
        handles = list(handles) + [_SpacerLegendHandle()] + list(extra_handles)
        labels = list(labels) + [""] + list(extra_labels or [""] * len(extra_handles))
        kwargs["labelspacing"] = 1.15
        kwargs["handler_map"] = {
            _SpacerLegendHandle: _HandlerSpacerLegend(),
            _PulseLegendHandle: _HandlerPulseLegend(),
            **(handler_map or {}),
        }
    ax.legend(handles, labels, **kwargs)


class _SpacerLegendHandle:
    """Invisible legend row that adds vertical space before the pulse box."""


class _HandlerSpacerLegend(HandlerBase):
    """Reserve a short blank row in the legend."""

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        # Tiny invisible patch so matplotlib keeps the row height.
        pad = Rectangle((xdescent, ydescent), width, 0.55 * height,
                        facecolor="none", edgecolor="none", transform=trans)
        return [pad]


class _PulseLegendHandle:
    """Legend proxy: white rectangle with gray outline and ``text`` inside."""

    def __init__(self, text: str) -> None:
        self.text = text


class _HandlerPulseLegend(HandlerBase):
    """Draw ``_PulseLegendHandle`` as a white box with gray outline containing its text."""

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        label = orig_handle.text
        box_w = max(width * 2.8, 0.55 * fontsize * max(len(label), 1))
        box_h = height * 2.2
        # Drop the box a bit below the default legend baseline for breathing room.
        x0, y0 = xdescent, ydescent - 1.1 * height
        rect = Rectangle((x0, y0), box_w, box_h, facecolor="white", edgecolor="0.45",
                         linewidth=0.9, transform=trans)
        txt = MplText(x0 + 0.5 * box_w, y0 + 0.5 * box_h, label,
                      fontsize=max(6.0, fontsize - 1), ha="center", va="center",
                      transform=trans, color="0.15")
        return [rect, txt]


def _finish(fig: plt.Figure, axes: np.ndarray, plot_cfg: Dict[str, Any], common: Dict[str, Any],
            out_dir: Path, stem: str) -> Path:
    """Date axis, optional figure title, and save."""
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    axes[-1].set_xlabel("Date (2026)")
    title = plot_cfg.get("title")
    if title:
        fig.suptitle(title)
    return _save(fig, out_dir / f"{stem}.{common.get('format', 'png')}", common)


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


def measured_relative_water(lab: Dict[str, Any], porosity: float) -> pd.DataFrame:
    """Measured relative water content: each side's shallow and deep VWC divided by porosity.

    Columns are ``{side}_shallow`` and ``{side}_deep``. This is the same ``s = theta / n``
    the model stores. Deep probes are NaN from ``DEEP_PROBE_CRASH`` onward.
    """
    df = pd.read_csv(REPO_ROOT / lab["vwc_file"], parse_dates=["datetime"]).set_index("datetime").sort_index()
    out = pd.DataFrame(index=df.index)
    for side in lab["sides"]:
        out[f"{side}_shallow"] = df[f"tree{lab['tree']}_{side}_shallow"] / porosity
        deep = df[f"tree{lab['tree']}_{side}_deep"] / porosity
        deep = deep.where(deep.index < DEEP_PROBE_CRASH)
        out[f"{side}_deep"] = deep
    return out


def _load_irrigation(lab: Dict[str, Any], start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    """Irrigation / salt pulse table clipped to ``[start, end]``."""
    irrigation = pd.read_csv(REPO_ROOT / lab["irrigation_file"], parse_dates=["timestamp"]).set_index("timestamp")
    return irrigation[(irrigation.index >= start) & (irrigation.index <= end)]


def _irrigation_days(irrigation: pd.DataFrame, lab: Dict[str, Any], side: str) -> pd.DatetimeIndex:
    """Calendar days when ``side`` received a water pulse (vol_liters > 0)."""
    vol = irrigation[f"tree{lab['tree']}_{side}_vol_liters"]
    return vol[vol > 0].index.normalize().unique()


def _saline_irrigation_days(irrigation: pd.DataFrame, lab: Dict[str, Any], side: str) -> set:
    """Calendar days when ``side`` received both water and salt (moles > 0)."""
    prefix = f"tree{lab['tree']}_{side}"
    both = irrigation[(irrigation[f"{prefix}_vol_liters"] > 0) & (irrigation[f"{prefix}_moles"] > 0)]
    return set(both.index.normalize())


def _annotate_pulse(ax: plt.Axes, t: pd.Timestamp, note: str) -> None:
    """Vertical dashed pulse marker; white label with gray outline to the right of the line."""
    ax.axvline(t, color="0.4", ls=":", lw=0.8, zorder=2)
    ax.annotate(
        note, xy=(t, 0.92), xycoords=("data", "axes fraction"),
        xytext=(8, 0), textcoords="offset points",
        ha="left", va="top", fontsize="x-small",
        bbox=dict(boxstyle="square,pad=0.25", facecolor="white", edgecolor="0.45", alpha=0.95),
        zorder=5, clip_on=False,
    )


def _pulse_legend_kwargs(plot_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Extra legend args when ``plots.*.pulse_legend`` is set; otherwise empty."""
    text = plot_cfg.get("pulse_legend")
    if not text:
        return {}
    return {
        "extra_handles": [_PulseLegendHandle(text)],
        "extra_labels": [""],
        "handler_map": {_PulseLegendHandle: _HandlerPulseLegend()},
    }


def _draw_irrigation_day(ax: plt.Axes, day: pd.Timestamp, ylim: tuple[float, float],
                         label: str = "IRRIGATION") -> None:
    """Blank the irrigation calendar day with a black-outlined band and vertical label."""
    x0 = mdates.date2num(day) - 0.5
    ymin, ymax = ylim
    ax.add_patch(Rectangle(
        (x0, ymin), 1.0, ymax - ymin, fill=False, edgecolor="black", linewidth=1.2, zorder=6,
    ))
    ax.text(
        x0 + 0.5, 0.5 * (ymin + ymax), label, rotation=90,
        ha="center", va="center", fontsize=7, fontweight="bold", zorder=7,
        linespacing=1.1,
    )


def plot_sap_flux_vs_uptake(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                            common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: model root water uptake (qs) per compartment against sap flux on a ground-area
    basis for the matching side, one panel per side. Dotted lines mark irrigation, with salt noted."""
    lab = cfg["lab"]
    sap = lab_inputs.sap_flux_ground_basis(lab)
    sap = sap[(sap.index >= df.index[0]) & (sap.index <= df.index[-1])]
    irrigation = _load_irrigation(lab, df.index[0], df.index[-1])
    linewidth = float(common.get("linewidth", 1.2))
    titles = plot_cfg.get("panel_titles") or lab["sides"]

    width, height = common.get("figsize", [9, 3.0])
    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])),
                             sharex=True, gridspec_kw={"hspace": 0.45})
    for i, (ax, side) in enumerate(zip(axes, lab["sides"])):
        ax.plot(df.index, _compartment(df, "qs", i), color=common.get("model_color"), lw=linewidth,
                label=f"Model root uptake $q_{{s,{i + 1}}}$")
        ax.plot(sap.index, sap[side], color=common.get("data_color"), lw=linewidth, marker=".", ms=3,
                label=f"Sap flux, {side} root")
        prefix = f"tree{lab['tree']}_{side}"
        for t, row in irrigation[irrigation[f"{prefix}_vol_liters"] > 0].iterrows():
            note = f"{row[f'{prefix}_vol_liters']} L"
            if row[f"{prefix}_moles"] > 0:
                note += f", {row[f'{prefix}_moles']} mol"
            _annotate_pulse(ax, t, note)
        ax.set_title(titles[i], fontsize="medium")
        ax.set_ylabel(plot_cfg.get("ylabel", "Flux ($\\mu$m s$^{-1}$)"))
        ylim = _panel_ylim(plot_cfg.get("ylim"), i)
        if ylim is not None:
            ax.set_ylim(ylim)
        _legend(ax, common, **_pulse_legend_kwargs(plot_cfg))

    return _finish(fig, axes, plot_cfg, common, out_dir, "sap_flux_vs_uptake")


def plot_daily_uptake_vs_sap_flux(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                                  common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: daily root uptake against daily sap flux (mm of water), one panel per side.

    Both fluxes are ground-area rates. Negative root uptake is left out of the daily total;
    it still appears on the sap-flux time series. Shallow and deep probe VWC depletion
    (from ``vwc_flux.csv``, ground-area um/s) are added as extra bars. Irrigation days omit
    bars and show an IRRIGATION band; deep-probe bars stop at ``DEEP_PROBE_CRASH``.
    """
    lab = cfg["lab"]
    sap = lab_inputs.sap_flux_ground_basis(lab)
    sap = sap[(sap.index >= df.index[0]) & (sap.index <= df.index[-1])]
    vwc_flux_path = REPO_ROOT / lab.get("vwc_flux_file", "data/soil_data/vwc_data/vwc_flux.csv")
    vwc_flux = pd.read_csv(vwc_flux_path, parse_dates=["datetime"]).set_index("datetime").sort_index()
    vwc_flux = vwc_flux[(vwc_flux.index >= df.index[0]) & (vwc_flux.index <= df.index[-1])]
    irrigation = _load_irrigation(lab, df.index[0], df.index[-1])
    titles = plot_cfg.get("panel_titles") or lab["sides"]
    width, height = common.get("figsize", [9, 3.0])
    data_color = common.get("data_color")
    ylim_cfg = plot_cfg.get("ylim")

    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])), sharex=True)
    for i, (ax, side) in enumerate(zip(np.atleast_1d(axes), lab["sides"])):
        rwu = daily_water_mm(pd.Series(_compartment(df, "qs", i), index=df.index).clip(lower=0))
        measured = daily_water_mm(sap[side])
        prefix = f"tree{lab['tree']}_{side}"
        shallow = daily_water_mm(vwc_flux[f"{prefix}_shallow_vol_flux_scaled_to_root"])
        deep_flux = vwc_flux[f"{prefix}_deep_vol_flux_scaled_to_root"].where(
            vwc_flux.index < DEEP_PROBE_CRASH
        )
        deep = daily_water_mm(deep_flux)
        # No deep daily bar on/after the crash calendar day.
        deep = deep[deep.index < DEEP_PROBE_CRASH.normalize()]
        irrig_days = _irrigation_days(irrigation, lab, side)
        saline_days = _saline_irrigation_days(irrigation, lab, side)
        rwu, measured, shallow, deep = rwu.copy(), measured.copy(), shallow.copy(), deep.copy()
        for series in (rwu, measured, shallow, deep):
            series.loc[series.index.isin(irrig_days)] = np.nan
        days = rwu.index.union(measured.index).union(shallow.index).union(deep.index)
        days = days.union(irrig_days).sort_values()
        x = mdates.date2num(days)
        bar_width = 0.2
        offsets = (-1.5 * bar_width, -0.5 * bar_width, 0.5 * bar_width, 1.5 * bar_width)
        ax.bar(x + offsets[0], rwu.reindex(days), width=bar_width, color=common.get("model_color"),
               label="Model root uptake")
        ax.bar(x + offsets[1], measured.reindex(days), width=bar_width, color=data_color,
               alpha=0.85, label="Sap flux")
        ax.bar(x + offsets[2], shallow.reindex(days), width=bar_width, color=data_color,
               alpha=0.55, hatch="//", label="Shallow VWC")
        ax.bar(x + offsets[3], deep.reindex(days), width=bar_width, color=data_color,
               alpha=0.55, hatch="..", label="Deep VWC")
        ylim = _panel_ylim(ylim_cfg, i)
        if ylim is None:
            ylim = (-1.0, 3.0)
        ax.set_ylim(ylim)
        for day in irrig_days:
            band = "IRRIGATION\nWITH SALINE WATER" if day in saline_days else "IRRIGATION"
            _draw_irrigation_day(ax, day, tuple(ylim), label=band)
        ax.set_title(titles[i], fontsize="medium")
        ax.set_ylabel(plot_cfg.get("ylabel", "Daily total (mm)"))
        _legend(ax, common)
    axes[-1].xaxis_date()
    return _finish(fig, axes, plot_cfg, common, out_dir, "daily_uptake_vs_sap_flux")


def plot_relative_water_content(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                                common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: modeled relative water content against shallow and deep probe values.

    One panel per side. Measured series are volumetric water content divided by porosity,
    the same relative water content ``s`` the model reports for that compartment.
    """
    sim_cfg, lab = cfg["simulation"], cfg["lab"]
    porosity = getattr(soil, sim_cfg["soil_texture"])().N
    obs = measured_relative_water(lab, porosity)
    obs = obs[(obs.index >= df.index[0]) & (obs.index <= df.index[-1])]
    titles = plot_cfg.get("panel_titles") or lab["sides"]
    linewidth = float(common.get("linewidth", 1.2))
    data_color = common.get("data_color")
    width, height = common.get("figsize", [9, 3.0])

    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])), sharex=True)
    for i, (ax, side) in enumerate(zip(np.atleast_1d(axes), lab["sides"])):
        ax.plot(df.index, _compartment(df, "s", i), color=common.get("model_color"), lw=linewidth,
                label="Model $s$")
        ax.plot(obs.index, obs[f"{side}_shallow"], color=data_color, lw=linewidth, ls=":",
                label="Shallow")
        ax.plot(obs.index, obs[f"{side}_deep"], color=data_color, lw=linewidth, ls="--",
                label="Deep")
        ax.set_title(titles[i], fontsize="medium")
        ax.set_ylabel(plot_cfg.get("ylabel", "Relative water content $s$ (-)"))
        ylim = _panel_ylim(plot_cfg.get("ylim"), i)
        if ylim is not None:
            ax.set_ylim(ylim)
        _legend(ax, common)
    return _finish(fig, axes, plot_cfg, common, out_dir, "relative_water_content")


def plot_soil_salt(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                   common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: modeled soil salt concentration against measured pore-water concentration.

    One panel per side. Measured values are the pore-water columns named in
    ``lab.conductivity_cols`` (mM). Dotted lines mark a salt addition.
    """
    lab = cfg["lab"]
    obs = pd.read_csv(REPO_ROOT / lab["conductivity_file"], parse_dates=["datetime"]).set_index("datetime")
    obs = obs[(obs.index >= df.index[0]) & (obs.index <= df.index[-1])]
    irrigation = _load_irrigation(lab, df.index[0], df.index[-1])
    titles = plot_cfg.get("panel_titles") or lab["sides"]
    linewidth = float(common.get("linewidth", 1.2))
    width, height = common.get("figsize", [9, 3.0])

    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])),
                             sharex=True, gridspec_kw={"hspace": 0.45})
    for i, (ax, side) in enumerate(zip(np.atleast_1d(axes), lab["sides"])):
        ax.plot(df.index, _compartment(df, "cs", i), color=common.get("model_color"), lw=linewidth,
                label="Model soil salt")
        ax.plot(obs.index, obs[lab["conductivity_cols"][side]], color=common.get("data_color"), lw=linewidth,
                label="Measured pore water")
        moles = irrigation[f"tree{lab['tree']}_{side}_moles"]
        for t, added in moles[moles > 0].items():
            _annotate_pulse(ax, t, f"{added} mol")
        ax.set_title(titles[i], fontsize="medium")
        ax.set_ylabel(plot_cfg.get("ylabel", "Salt concentration (mM)"))
        ylim = _panel_ylim(plot_cfg.get("ylim"), i)
        if ylim is not None:
            ax.set_ylim(ylim)
        _legend(ax, common, **_pulse_legend_kwargs(plot_cfg))
    return _finish(fig, axes, plot_cfg, common, out_dir, "soil_salt")


def plot_matric_potential(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                          common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: modeled matric potential against the Board 9 matric probes.

    Modeled matric potential is the soil retention curve at the modeled relative
    water content. ``lab.matric_cols`` names the measured column for each side.
    """
    sim_cfg, lab = cfg["simulation"], cfg["lab"]
    texture = getattr(soil, sim_cfg["soil_texture"])()
    formulation = sim_cfg.get("soil_formulation")
    if formulation is not None and hasattr(texture, "FORMULATION"):
        texture.FORMULATION = str(formulation).strip().lower().replace("-", "_")
    obs = pd.read_csv(REPO_ROOT / lab["conductivity_file"], parse_dates=["datetime"]).set_index("datetime")
    obs = obs[(obs.index >= df.index[0]) & (obs.index <= df.index[-1])]
    titles = plot_cfg.get("panel_titles") or lab["sides"]
    linewidth = float(common.get("linewidth", 1.2))
    width, height = common.get("figsize", [9, 3.0])

    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])), sharex=True)
    for i, (ax, side) in enumerate(zip(np.atleast_1d(axes), lab["sides"])):
        modeled = np.asarray(texture.psi_s(_compartment(df, "s", i)), dtype=float)
        ax.plot(df.index, modeled, color=common.get("model_color"), lw=linewidth,
                label="Model matric potential")
        ax.plot(obs.index, obs[lab["matric_cols"][side]], color=common.get("data_color"), lw=linewidth,
                label="Measured matric potential")
        ax.set_title(titles[i], fontsize="medium")
        ax.set_ylabel(plot_cfg.get("ylabel", "Matric potential (MPa)"))
        ylim = _panel_ylim(plot_cfg.get("ylim"), i)
        if ylim is not None:
            ax.set_ylim(ylim)
        _legend(ax, common)
    return _finish(fig, axes, plot_cfg, common, out_dir, "matric_potential")


# Config ``plots`` key -> plotting function.
PLOT_FUNCTIONS: Dict[str, Callable[..., Path]] = {
    "sap_flux_vs_uptake": plot_sap_flux_vs_uptake,
    "daily_uptake_vs_sap_flux": plot_daily_uptake_vs_sap_flux,
    "relative_water_content": plot_relative_water_content,
    "soil_salt": plot_soil_salt,
    "matric_potential": plot_matric_potential,
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
