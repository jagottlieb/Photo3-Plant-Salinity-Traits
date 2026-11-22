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
    python lab_experiment_simulation/post_processing.py lab_experiment_simulation/salt_treatment/config.yaml --only salt_concentration_attribution
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
import matplotlib.ticker as mticker
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import yaml

# Deep TEROS probes begin failing ~16:00 on 8/17 and are gone by 17:00 (both sides).
DEEP_PROBE_CRASH = pd.Timestamp("2026-08-17 16:00:00")

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "sensitivity_analysis"))
sys.path.append(str(REPO_ROOT))

from post_processing import _panel_ylim, set_plot_style  # sensitivity_analysis/post_processing.py, first on sys.path

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
        # Keep outside legends inside the saved tight bbox (pulse box + long labels).
        "borderaxespad": float(legend_cfg.get("borderaxespad", 0.5)),
    }
    if extra_handles:
        handles = list(handles) + list(extra_handles)
        labels = list(labels) + list(extra_labels or [""] * len(extra_handles))
        kwargs["handler_map"] = {
            _PulseLegendHandle: _HandlerPulseLegend(),
            **(handler_map or {}),
        }
    ax.legend(handles, labels, **kwargs)


class _PulseLegendHandle:
    """Legend proxy: white rectangle with gray outline and ``text`` inside."""

    def __init__(self, text: str) -> None:
        self.text = text


class _HandlerPulseLegend(HandlerBase):
    """Draw ``_PulseLegendHandle`` as a compact white/gray box (label text sits beside it).

    Keeping the box within the normal legend handle width so ``bbox_inches='tight'``
    does not clip it; the descriptive string is the legend label, not text inside the box.
    """

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        box_w = max(width, 1.6 * height)
        box_h = 1.35 * height
        x0 = xdescent + 0.5 * (width - box_w)
        y0 = ydescent + 0.5 * (height - box_h)
        rect = Rectangle((x0, y0), box_w, box_h, facecolor="white", edgecolor="0.45",
                         linewidth=0.9, transform=trans)
        return [rect]


def _apply_date_axis(ax: plt.Axes, common: Dict[str, Any]) -> None:
    """Daily vertical gridlines and labels every N days; no x tick marks.

    ``plots.common.date_label_interval_days`` (default 2) sets the label spacing.
    ``plots.common.date_label_anchor`` (e.g. ``2026-08-12``) pins the labeled
    days and, when set, also sets the left edge of the x-axis.
    ``plots.common.xlim_pad_before_days`` (default 0) extends that left edge
    earlier so day-centered irrigation bands (half-day before midnight) are not
    clipped — e.g. ``0.6`` for the salt-treatment daily bar plots.

    Vertical day lines are drawn with ``axvline`` because seaborn's minor x-grid
    does not reliably render under the Agg + ``bbox_inches='tight'`` save path.
    """
    interval = int(common.get("date_label_interval_days", 2))
    pad_before = float(common.get("xlim_pad_before_days", 0.0))
    _, x1 = ax.get_xlim()
    if common.get("date_label_anchor") is not None:
        anchor = pd.Timestamp(common["date_label_anchor"]).normalize()
        ax.set_xlim(left=mdates.date2num(anchor) - pad_before, right=x1)
    else:
        anchor = pd.Timestamp(mdates.num2date(ax.get_xlim()[0])).tz_localize(None).normalize()
        if pad_before:
            ax.set_xlim(left=mdates.date2num(anchor) - pad_before, right=x1)
    x0, x1 = ax.get_xlim()
    d0 = pd.Timestamp(mdates.num2date(x0)).tz_localize(None).normalize()
    d1 = pd.Timestamp(mdates.num2date(x1)).tz_localize(None).normalize()
    offset = (d0 - anchor).days % interval
    first = d0 if offset == 0 else d0 + pd.Timedelta(days=interval - offset)
    majors = list(pd.date_range(first, d1, freq=f"{interval}D")) if first <= d1 else []
    ax.xaxis.set_major_locator(mticker.FixedLocator(mdates.date2num(majors)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.tick_params(axis="x", which="both", length=0, bottom=False)
    grid_color = ax.xaxis.get_gridlines()[0].get_color() if ax.xaxis.get_gridlines() else "0.8"
    ax.grid(False, axis="x")
    ax.set_axisbelow(True)
    for day in pd.date_range(d0, d1, freq="D"):
        ax.axvline(day, color=grid_color, lw=0.8, zorder=1)


def _finish(fig: plt.Figure, axes: np.ndarray, plot_cfg: Dict[str, Any], common: Dict[str, Any],
            out_dir: Path, stem: str) -> Path:
    """Date axis, optional figure title, and save."""
    axes_arr = np.atleast_1d(axes)
    for ax in axes_arr:
        _apply_date_axis(ax, common)
    axes_arr[-1].set_xlabel("Date (2026)")
    title = plot_cfg.get("title")
    if title:
        fig.suptitle(title)
    # Extra pad so outside legends (pulse box + labels) are not clipped by tight bbox.
    pad = float(common.get("savefig_pad_inches", 0.25))
    out_path = out_dir / f"{stem}.{common.get('format', 'png')}"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=common.get("dpi", 300), bbox_inches="tight", pad_inches=pad)
    plt.close(fig)
    return out_path


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


def daily_vwc_depletion_mm(theta: pd.Series, volume_mL: float, ground_area_cm2: float) -> pd.Series:
    """Daily soil-water loss (mm) from midnight-to-midnight VWC change.

    For calendar day ``D``: ``-(θ(D+1 00:00) - θ(D 00:00)) * V / A * 1000``.
    Depletion is positive. Uses the whole pot volume so the estimate is on the same
    ground-area basis as integrated root uptake / sap flux (mm). Requires exact
    midnight samples in ``theta``.
    """
    at_midnight = theta.dropna().sort_index()
    at_midnight = at_midnight.loc[at_midnight.index == at_midnight.index.normalize()]
    at_midnight = at_midnight[~at_midnight.index.duplicated(keep="first")]
    if len(at_midnight) < 2:
        return pd.Series(dtype=float)
    vol_m3 = float(volume_mL) * 1e-6
    area_m2 = float(ground_area_cm2) * 1e-4
    scale = vol_m3 / area_m2 * 1000.0
    out: Dict[pd.Timestamp, float] = {}
    for t0, v0 in at_midnight.items():
        t1 = t0 + pd.Timedelta(days=1)
        if t1 not in at_midnight.index:
            continue
        out[t0] = -float(at_midnight[t1] - v0) * scale
    return pd.Series(out, dtype=float)


def daily_soil_evaporation_mm(s: pd.Series, evmax_mm_d: float, sh: float) -> pd.Series:
    """Fallback daily soil evaporation (mm) reconstructed from ``s(t)``.

    Used only when the timeseries lacks ``soil_evap`` (older batches). Prefer
    integrating the model ``soil_evap`` column (um/s) with ``daily_water_mm``.
    Uses ``EVMAX * (s - SH) / (1 - SH)`` mm/day when ``s > SH``, else 0.
    """
    s = s.dropna().sort_index()
    if len(s) < 2:
        return pd.Series(dtype=float)
    dt = s.index.to_series().diff().dt.total_seconds()
    dt.iloc[0] = dt.iloc[1]
    s_vals = s.to_numpy(dtype=float)
    evap_mm_d = np.where(s_vals > sh, evmax_mm_d * (s_vals - sh) / (1.0 - sh), 0.0)
    depth_mm = evap_mm_d * dt.to_numpy() / 86400.0
    day = (s.index - pd.Timedelta(seconds=1)).normalize()
    return pd.Series(depth_mm, index=day).groupby(level=0).sum()


def daily_model_soil_evaporation_mm(df: pd.DataFrame, i: int, texture: Any) -> pd.Series:
    """Daily soil evaporation (mm) for compartment ``i`` from the model timeseries.

    Prefers ``soil_evap`` (um/s, written by ``SoilMultiple``). Falls back to
    reconstructing from ``s`` with the texture ``EVMAX`` / ``SH`` for older batches.
    """
    if "soil_evap" in df.columns:
        return daily_water_mm(pd.Series(_compartment(df, "soil_evap", i), index=df.index))
    evmax = float(getattr(texture, "EVMAX", soil.SoilMultiple.EVMAX))
    sh = float(texture.SH)
    return daily_soil_evaporation_mm(pd.Series(_compartment(df, "s", i), index=df.index), evmax, sh)


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


def _annotate_pulse(ax: plt.Axes, t: pd.Timestamp, note: str, placement: str = "in_pane",
                    y: float = 0.92) -> None:
    """Vertical dashed pulse marker with a white / gray-outline label.

    ``placement`` is ``in_pane`` (to the right of the line) or ``above`` (above the axes).
    For ``in_pane``, ``y`` is the label position in axes fraction (1 = top).
    """
    ax.axvline(t, color="0.4", ls=":", lw=0.8, zorder=2)
    bbox = dict(boxstyle="square,pad=0.25", facecolor="white", edgecolor="0.45", alpha=0.95)
    if placement == "above":
        ax.text(
            t, 1.08, note, transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize="x-small", bbox=bbox,
            zorder=5, clip_on=False,
        )
    else:
        ax.annotate(
            note, xy=(t, y), xycoords=("data", "axes fraction"),
            xytext=(8, 0), textcoords="offset points",
            ha="left", va="top", fontsize="x-small", bbox=bbox,
            zorder=5, clip_on=False,
        )


def _pulse_legend_kwargs(text: Optional[str]) -> Dict[str, Any]:
    """Extra legend args for a pulse-annotation box; empty if ``text`` is falsy."""
    if not text:
        return {}
    return {
        "extra_handles": [_PulseLegendHandle(text)],
        "extra_labels": [text],
        "handler_map": {_PulseLegendHandle: _HandlerPulseLegend()},
    }


def _pulse_legend_text(has_volume: bool, has_salt: bool) -> Optional[str]:
    """Legend label for pulse boxes on a panel, or None if that panel has no pulses."""
    if has_volume and has_salt:
        return "volume, salt mass added"
    if has_salt:
        return "salt mass added"
    if has_volume:
        return "volume added"
    return None


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


def _mean_rpd(model: pd.Series, observed: pd.Series) -> float:
    """Mean relative percent difference: 100 * mean(|m-o| / ((|m|+|o|)/2))."""
    denom = 0.5 * (model.abs() + observed.abs())
    valid = denom > 1e-12
    if not valid.any():
        return float("nan")
    return float(100.0 * np.mean((model - observed).abs()[valid] / denom[valid]))


def _write_metrics_csv(path: Path, rows: List[Dict[str, Any]]) -> Path:
    """Write ``rows`` as a CSV next to the figures; floats rounded to 3 decimals."""
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows)
    float_cols = table.select_dtypes(include=["floating"]).columns
    if len(float_cols):
        table[float_cols] = table[float_cols].round(3)
    table.to_csv(path, index=False)
    print(f"  Wrote {path.name}")
    return path


def plot_sap_flux_vs_uptake(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                            common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: model root water uptake (qs) per compartment against sap flux on a ground-area
    basis for the matching side, one panel per side. Dotted lines mark irrigation, with salt noted.

    When ``show_error_metrics`` is true, writes ``sap_flux_vs_uptake_metrics.csv`` (RMSE, R²,
    peak RPD from daily maxima, overall RPD) with no on-figure annotations.
    """
    lab = cfg["lab"]
    sap = lab_inputs.sap_flux_ground_basis(lab)
    sap = sap[(sap.index >= df.index[0]) & (sap.index <= df.index[-1])]
    irrigation = _load_irrigation(lab, df.index[0], df.index[-1])
    linewidth = float(common.get("linewidth", 1.2))
    titles = plot_cfg.get("panel_titles") or lab["sides"]
    export_metrics = bool(plot_cfg.get("show_error_metrics", False))
    metric_rows: List[Dict[str, Any]] = []

    width, height = common.get("figsize", [9, 3.0])
    placement = str(plot_cfg.get("pulse_placement", "in_pane"))
    pulse_y = float(plot_cfg.get("pulse_y", 0.92))
    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])),
                             sharex=True, gridspec_kw={"hspace": 0.45})
    for i, (ax, side) in enumerate(zip(axes, lab["sides"])):
        qs = pd.Series(_compartment(df, "qs", i), index=df.index)
        ax.plot(df.index, qs, color=common.get("model_color"), lw=linewidth,
                label=f"Model root uptake $q_{{s,{i + 1}}}$")
        ax.plot(sap.index, sap[side], color=common.get("data_color"), lw=linewidth, marker=".", ms=3,
                label=f"Sap flux, {side} root")
        if export_metrics:
            paired = pd.concat([qs.rename("qs"), sap[side].rename("sap")], axis=1, join="inner").dropna()
            if not paired.empty:
                err = paired["qs"] - paired["sap"]
                rmse = float(np.sqrt(np.mean(err**2)))
                ss_tot = float(np.sum((paired["sap"] - paired["sap"].mean())**2))
                r2 = float(1.0 - np.sum(err**2) / ss_tot) if ss_tot > 0 else float("nan")
                peaks = paired.groupby(paired.index.normalize()).max()
                metric_rows.append({
                    "side": side,
                    "panel_title": titles[i],
                    "rmse_um_s": rmse,
                    "r2": r2,
                    "peak_rpd_pct": _mean_rpd(peaks["qs"], peaks["sap"]),
                    "overall_rpd_pct": _mean_rpd(paired["qs"], paired["sap"]),
                    "n_paired": int(len(paired)),
                })
        prefix = f"tree{lab['tree']}_{side}"
        has_volume = has_salt = False
        for t, row in irrigation[irrigation[f"{prefix}_vol_liters"] > 0].iterrows():
            note = f"{row[f'{prefix}_vol_liters']} L"
            has_volume = True
            if row[f"{prefix}_moles"] > 0:
                note += f", {row[f'{prefix}_moles']} mol"
                has_salt = True
            _annotate_pulse(ax, t, note, placement=placement, y=pulse_y)
        ax.set_title(titles[i], fontsize="medium")
        ax.set_ylabel(plot_cfg.get("ylabel", "Flux ($\\mu$m s$^{-1}$)"))
        ylim = _panel_ylim(plot_cfg.get("ylim"), i)
        if ylim is not None:
            ax.set_ylim(ylim)
        _legend(ax, common, **_pulse_legend_kwargs(_pulse_legend_text(has_volume, has_salt)))

    if export_metrics and metric_rows:
        _write_metrics_csv(out_dir / "sap_flux_vs_uptake_metrics.csv", metric_rows)
    return _finish(fig, axes, plot_cfg, common, out_dir, "sap_flux_vs_uptake")


def plot_daily_uptake_vs_sap_flux(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                                  common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: daily root uptake against daily sap flux (mm of water), one panel per side.

    Model root uptake is optionally stacked with model soil evaporation (gray) when
    ``include_soil_evaporation`` is true (default). Evaporation comes from the model
    ``soil_evap`` timeseries (um/s); older batches without that column fall back to a
    reconstruction from ``s(t)``. Sap compares to root uptake; probe bars
    (midnight-to-midnight VWC → whole-pot mm) compare to total model soil loss
    (RWU, or RWU+Ev when evaporation is included). Irrigation days omit bars; deep
    probes stop before ``DEEP_PROBE_CRASH``.

    When ``show_error_metrics`` is true, writes daily / summary CSVs with evaporation and
    total model soil flux, and differences vs sap (RWU) or VWC (RWU[+Ev]).
    """
    lab = cfg["lab"]
    sim_cfg = cfg["simulation"]
    day_lo = pd.Timestamp(sim_cfg["start_date"]).normalize()
    day_hi = pd.Timestamp(sim_cfg["end_date"]).normalize()
    sap = lab_inputs.sap_flux_ground_basis(lab)
    sap = sap[(sap.index >= df.index[0]) & (sap.index <= df.index[-1])]
    vwc_path = REPO_ROOT / lab["vwc_file"]
    vwc = pd.read_csv(vwc_path, parse_dates=["datetime"]).set_index("datetime").sort_index()
    # Calendar midnights for [start_date, end_date+1] so each comparison day has a closing sample.
    vwc = vwc[(vwc.index >= day_lo) & (vwc.index <= day_hi + pd.Timedelta(days=1))]
    irrigation = _load_irrigation(lab, df.index[0], df.index[-1])
    titles = plot_cfg.get("panel_titles") or lab["sides"]
    export_metrics = bool(plot_cfg.get("show_error_metrics", False))
    include_evap = bool(plot_cfg.get("include_soil_evaporation", True))
    width, height = common.get("figsize", [9, 3.0])
    data_color = common.get("data_color")
    evap_color = common.get("evap_color", "#969696")
    ylim_cfg = plot_cfg.get("ylim")
    ground_area_cm2 = float(lab["ground_area_cm2"])
    texture = getattr(soil, sim_cfg["soil_texture"])()
    daily_rows: List[Dict[str, Any]] = []
    summary_rows: List[Dict[str, Any]] = []

    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])), sharex=True)
    for i, (ax, side) in enumerate(zip(np.atleast_1d(axes), lab["sides"])):
        qs = pd.Series(_compartment(df, "qs", i), index=df.index)
        rwu = daily_water_mm(qs.clip(lower=0))
        if include_evap:
            evap = daily_model_soil_evaporation_mm(df, i, texture)
        else:
            evap = pd.Series(0.0, index=rwu.index, dtype=float)
        measured = daily_water_mm(sap[side])
        prefix = f"tree{lab['tree']}_{side}"
        volume_mL = float(lab["container_volume_mL"][side])
        shallow = daily_vwc_depletion_mm(vwc[f"{prefix}_shallow"], volume_mL, ground_area_cm2)
        deep_theta = vwc[f"{prefix}_deep"].where(vwc.index <= DEEP_PROBE_CRASH.normalize())
        deep = daily_vwc_depletion_mm(deep_theta, volume_mL, ground_area_cm2)
        # Drop the crash calendar day and after (incomplete / failed deep record).
        deep = deep[deep.index < DEEP_PROBE_CRASH.normalize()]
        irrig_days = _irrigation_days(irrigation, lab, side)
        saline_days = _saline_irrigation_days(irrigation, lab, side)
        rwu, evap, measured, shallow, deep = (
            rwu.copy(), evap.copy(), measured.copy(), shallow.copy(), deep.copy()
        )
        for series in (rwu, evap, measured, shallow, deep):
            series.loc[series.index.isin(irrig_days)] = np.nan
        days = rwu.index.union(evap.index).union(measured.index).union(shallow.index).union(deep.index)
        days = days.union(irrig_days).sort_values()
        days = days[(days >= day_lo) & (days <= day_hi)]
        rwu_d = rwu.reindex(days)
        evap_d = evap.reindex(days) if include_evap else pd.Series(np.nan, index=days, dtype=float)
        total_d = rwu_d + evap_d if include_evap else rwu_d
        sap_d = measured.reindex(days)
        shallow_d = shallow.reindex(days)
        deep_d = deep.reindex(days)
        if export_metrics:
            for day in days:
                model_mm = float(rwu_d.loc[day]) if pd.notna(rwu_d.loc[day]) else np.nan
                evap_mm = float(evap_d.loc[day]) if include_evap and pd.notna(evap_d.loc[day]) else (
                    0.0 if not include_evap and pd.notna(rwu_d.loc[day]) else np.nan
                )
                total_mm = float(total_d.loc[day]) if pd.notna(total_d.loc[day]) else np.nan
                sap_mm = float(sap_d.loc[day]) if pd.notna(sap_d.loc[day]) else np.nan
                shallow_mm = float(shallow_d.loc[day]) if pd.notna(shallow_d.loc[day]) else np.nan
                deep_mm = float(deep_d.loc[day]) if pd.notna(deep_d.loc[day]) else np.nan
                daily_rows.append({
                    "date": pd.Timestamp(day).strftime("%Y-%m-%d"),
                    "side": side,
                    "panel_title": titles[i],
                    "is_irrigation_day": bool(day in irrig_days),
                    "model_rwu_mm": model_mm,
                    "model_evap_mm": evap_mm,
                    "model_total_soil_mm": total_mm,
                    "sap_mm": sap_mm,
                    "shallow_vwc_mm": shallow_mm,
                    "deep_vwc_mm": deep_mm,
                    "model_rwu_minus_sap_mm": model_mm - sap_mm,
                    "model_total_minus_shallow_mm": total_mm - shallow_mm,
                    "model_total_minus_deep_mm": total_mm - deep_mm,
                })
            diff_sap = rwu_d - sap_d
            diff_shallow = total_d - shallow_d
            diff_deep = total_d - deep_d
            summary_rows.append({
                "side": side,
                "panel_title": titles[i],
                "n_days_rwu_vs_sap": int(diff_sap.notna().sum()),
                "n_days_total_vs_shallow": int(diff_shallow.notna().sum()),
                "n_days_total_vs_deep": int(diff_deep.notna().sum()),
                "mean_model_rwu_minus_sap_mm": float(diff_sap.mean()),
                "mean_model_total_minus_shallow_mm": float(diff_shallow.mean()),
                "mean_model_total_minus_deep_mm": float(diff_deep.mean()),
                "total_model_rwu_minus_sap_mm": float(diff_sap.sum()),
                "total_model_total_minus_shallow_mm": float(diff_shallow.sum()),
                "total_model_total_minus_deep_mm": float(diff_deep.sum()),
                "total_model_rwu_mm": float(rwu_d.sum()),
                "total_model_evap_mm": float(evap_d.fillna(0).sum()) if include_evap else 0.0,
                "total_model_total_soil_mm": float(total_d.sum()),
                "total_sap_mm": float(sap_d.sum()),
                "total_shallow_vwc_mm": float(shallow_d.sum()),
                "total_deep_vwc_mm": float(deep_d.sum()),
            })
        x = mdates.date2num(days)
        bar_width = 0.2
        offsets = (-1.5 * bar_width, -0.5 * bar_width, 0.5 * bar_width, 1.5 * bar_width)
        ax.bar(x + offsets[0], rwu_d, width=bar_width, color=common.get("model_color"),
               label="Model root uptake")
        if include_evap:
            ax.bar(x + offsets[0], evap_d, width=bar_width, bottom=rwu_d, color=evap_color,
                   label="Model soil evaporation")
        ax.bar(x + offsets[1], sap_d, width=bar_width, color=data_color,
               alpha=0.85, label="Sap flux")
        ax.bar(x + offsets[2], shallow_d, width=bar_width, color=data_color,
               alpha=0.55, hatch="//", label="Shallow VWC")
        ax.bar(x + offsets[3], deep_d, width=bar_width, color=data_color,
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
    if export_metrics and daily_rows:
        _write_metrics_csv(out_dir / "daily_uptake_vs_sap_flux_daily.csv", daily_rows)
        _write_metrics_csv(out_dir / "daily_uptake_vs_sap_flux_summary.csv", summary_rows)
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
    pulse_y = float(plot_cfg.get("pulse_y", 0.22))

    fig, axes = plt.subplots(len(lab["sides"]), 1, figsize=(width, height * len(lab["sides"])),
                             sharex=True, gridspec_kw={"hspace": 0.45})
    for i, (ax, side) in enumerate(zip(np.atleast_1d(axes), lab["sides"])):
        ax.plot(df.index, _compartment(df, "cs", i), color=common.get("model_color"), lw=linewidth,
                label="Model soil salt")
        ax.plot(obs.index, obs[lab["conductivity_cols"][side]], color=common.get("data_color"), lw=linewidth,
                label="Measured pore water")
        prefix = f"tree{lab['tree']}_{side}"
        vol = irrigation[f"{prefix}_vol_liters"]
        moles = irrigation[f"{prefix}_moles"]
        has_volume = has_salt = False
        for t in irrigation.index:
            v, m = float(vol.loc[t]), float(moles.loc[t])
            if v <= 0 and m <= 0:
                continue
            parts = []
            if v > 0:
                parts.append(f"{v} L")
                has_volume = True
            if m > 0:
                parts.append(f"{m} mol")
                has_salt = True
            _annotate_pulse(ax, t, ", ".join(parts), y=pulse_y)
        ax.set_title(titles[i], fontsize="medium")
        ax.set_ylabel(plot_cfg.get("ylabel", "Salt concentration (mM)"))
        ylim = _panel_ylim(plot_cfg.get("ylim"), i)
        if ylim is not None:
            ax.set_ylim(ylim)
        _legend(ax, common, **_pulse_legend_kwargs(_pulse_legend_text(has_volume, has_salt)))
    return _finish(fig, axes, plot_cfg, common, out_dir, "soil_salt")


def plot_salt_concentrations(df: pd.DataFrame, cfg: Dict[str, Any], plot_cfg: Dict[str, Any],
                             common: Dict[str, Any], out_dir: Path) -> Path:
    """Figure type: stem/leaf storage and left/right soil salt concentrations (mM).

    Storage concentrations are ``c_stem`` and ``c_leaf``. Soil concentrations are
    compartment ``cs`` (left dashed black, right solid black). Salt-addition
    pulses from the irrigation sheet are annotated.
    """
    lab = cfg["lab"]
    irrigation = _load_irrigation(lab, df.index[0], df.index[-1])
    linewidth = float(common.get("linewidth", 1.2))
    width, height = common.get("figsize", [9, 3.0])
    model_color = common.get("model_color", "#08519c")
    leaf_color = common.get("data_color", "#6baed6")
    pulse_y = float(plot_cfg.get("pulse_y", 0.92))
    title = plot_cfg.get("title")

    fig, ax = plt.subplots(figsize=(width, height))
    ax.plot(df.index, df["c_stem"], color=model_color, lw=linewidth,
            label=r"$c_{stem}$")
    ax.plot(df.index, df["c_leaf"], color=leaf_color, lw=linewidth,
            label=r"$c_{leaf}$")
    ax.plot(df.index, _compartment(df, "cs", 0), color="k", lw=linewidth, ls="--",
            label=r"Left soil $c_s$")
    ax.plot(df.index, _compartment(df, "cs", 1), color="k", lw=linewidth, ls="-",
            label=r"Right soil $c_s$")

    has_volume = has_salt = False
    for side in lab["sides"]:
        prefix = f"tree{lab['tree']}_{side}"
        for t, row in irrigation.iterrows():
            v, m = float(row[f"{prefix}_vol_liters"]), float(row[f"{prefix}_moles"])
            if m <= 0:
                continue
            note = f"{v} L, {m} mol" if v > 0 else f"{m} mol"
            has_salt = True
            if v > 0:
                has_volume = True
            _annotate_pulse(ax, t, note, y=pulse_y)

    ax.set_ylabel(plot_cfg.get("ylabel", "Salt concentration (mM)"))
    ylim = plot_cfg.get("ylim")
    if ylim is not None:
        ax.set_ylim(ylim)
    if title:
        fig.suptitle(title)
    _legend(ax, common, **_pulse_legend_kwargs(_pulse_legend_text(has_volume, has_salt)))
    return _finish(fig, np.array([ax]), plot_cfg, common, out_dir, "salt_concentrations")


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


def export_salt_concentration_attribution(
    df: pd.DataFrame,
    cfg: Dict[str, Any],
    export_cfg: Dict[str, Any],
    out_dir: Path,
) -> Path:
    """Write an Excel workbook of soil-salt concentration attribution by side.

    One sheet per root compartment. Columns ``s``, ``cs``, and ``MS`` are values
    from the model timeseries. Step changes and concentration terms are Excel
    formulas so the workbook shows the full calculation:

    - ``ds``, ``dMS`` from consecutive rows
    - ``dcs_actual = cs_t - cs_{t-1}``
    - ``dcs_uptake = dMS * cs_prev / MS_prev`` (salt-mass term)
    - ``dcs_s = -cs_prev * ds / s_prev`` (moisture term)
    - ``dcs_attributed = dcs_uptake + dcs_s`` (should match ``dcs_actual``)
    - ``dcs_residual = dcs_actual - dcs_attributed`` (finite-step truncation)
    - ``pct_abs_uptake_vs_abs_s = 100 * |dcs_uptake| / |dcs_s|``
    - ``is_pulse`` when ``ds`` or ``dMS`` exceeds ``pulse_threshold``

    A ``Summary`` sheet aggregates non-pulse steps with ``SUMIF``. A ``Notes``
    sheet documents units and formulas.
    """
    from openpyxl import Workbook
    from openpyxl.styles import Font
    from openpyxl.utils import get_column_letter

    lab = cfg["lab"]
    sides = list(lab["sides"])
    titles = list(export_cfg.get("panel_titles") or [])
    if len(titles) < len(sides):
        default_titles = {
            "left": "Left (Salt treated)",
            "right": "Right (Baseline)",
        }
        titles = [default_titles.get(side, side) for side in sides]
    pulse_threshold = float(export_cfg.get("pulse_threshold", 0.02))
    area_m2 = float(lab["ground_area_cm2"]) * 1e-4

    headers = [
        "datetime",
        "s",
        "cs_mol_m3",
        "MS_mol_m2",
        "ds",
        "dMS_mol_m2",
        "dMS_mol",
        "dcs_actual_mol_m3",
        "dcs_uptake_mol_m3",
        "dcs_s_mol_m3",
        "dcs_attributed_mol_m3",
        "dcs_residual_mol_m3",
        "abs_dcs_uptake",
        "abs_dcs_s",
        "pct_abs_uptake_vs_abs_s",
        "is_pulse",
    ]

    wb = Workbook()
    # Notes first for context when opening the file.
    ws_notes = wb.active
    ws_notes.title = "Notes"
    notes = [
        "Soil salt concentration attribution (model timeseries, burn-in excluded).",
        "",
        "Identity: cs = MS / (s * N * ZR)   [mol/m3]",
        "MS is salt mass per unit ground area [mol/m2].",
        f"Ground area A = {lab['ground_area_cm2']} cm2 = {area_m2} m2; dMS_mol = dMS_mol_m2 * A.",
        "",
        "Total differential (N*ZR fixed):",
        "  dcs = (partial cs / partial MS) dMS + (partial cs / partial s) ds",
        "      = dMS * (cs / MS)  -  cs * (ds / s)",
        "So at each timestep the two attributed terms should reconstruct dcs:",
        "  dcs_actual = cs_t - cs_{t-1}",
        "  dcs_uptake + dcs_s = dcs_attributed  ~  dcs_actual",
        "  dcs_residual = dcs_actual - dcs_attributed   (usually tiny; finite-step error)",
        "",
        "At each timestep after the first (using previous-row s, cs, MS):",
        "  ds = s_t - s_{t-1}",
        "  dMS_mol_m2 = MS_t - MS_{t-1}",
        "  dcs_actual_mol_m3 = cs_t - cs_{t-1}",
        "  dcs_uptake_mol_m3 = dMS_mol_m2 * cs_{t-1} / MS_{t-1}",
        "  dcs_s_mol_m3 = -cs_{t-1} * ds / s_{t-1}",
        "  dcs_attributed_mol_m3 = dcs_uptake + dcs_s",
        "  dcs_residual_mol_m3 = dcs_actual - dcs_attributed",
        "  pct_abs_uptake_vs_abs_s = 100 * |dcs_uptake| / |dcs_s|",
        f"  is_pulse = OR(ds > {pulse_threshold}, dMS_mol_m2 > {pulse_threshold})",
        "",
        "Negative dMS (and dcs_uptake) means salt mass left the soil (plant uptake when deplete_soil_salt is on).",
        "Positive dcs_s means concentration rose because soil moisture fell.",
        "Irrigation / salt pulses are flagged by is_pulse; Summary uses non-pulse rows only.",
        "",
        "Side sheets use Excel formulas for all derived columns so values recalculate if inputs change.",
    ]
    for i, line in enumerate(notes, start=1):
        ws_notes.cell(i, 1, line)
    ws_notes.column_dimensions["A"].width = 110

    def _write_side_sheet(title: str, side_index: int) -> int:
        ws = wb.create_sheet(title[:31])  # Excel sheet-name limit
        for col, name in enumerate(headers, start=1):
            cell = ws.cell(1, col, name)
            cell.font = Font(bold=True)

        s_vals = _compartment(df, "s", side_index)
        cs_vals = _compartment(df, "cs", side_index)
        ms_vals = _compartment(df, "MS", side_index)
        times = list(df.index)

        for r, (t, s, cs, ms) in enumerate(zip(times, s_vals, cs_vals, ms_vals), start=2):
            ws.cell(r, 1, t.to_pydatetime() if hasattr(t, "to_pydatetime") else t)
            ws.cell(r, 2, float(s))
            ws.cell(r, 3, float(cs))
            ws.cell(r, 4, float(ms))
            if r == 2:
                continue
            # Excel formulas (row r uses previous row r-1 as start-of-step state).
            ws.cell(r, 5, f"=B{r}-B{r-1}")  # ds
            ws.cell(r, 6, f"=D{r}-D{r-1}")  # dMS_mol_m2
            ws.cell(r, 7, f"=F{r}*{area_m2}")  # dMS_mol
            ws.cell(r, 8, f"=C{r}-C{r-1}")  # dcs_actual
            ws.cell(r, 9, f"=F{r}*C{r-1}/D{r-1}")  # dcs_uptake
            ws.cell(r, 10, f"=-C{r-1}*E{r}/B{r-1}")  # dcs_s
            ws.cell(r, 11, f"=I{r}+J{r}")  # dcs_attributed
            ws.cell(r, 12, f"=H{r}-K{r}")  # dcs_residual
            ws.cell(r, 13, f"=ABS(I{r})")  # abs_dcs_uptake
            ws.cell(r, 14, f"=ABS(J{r})")  # abs_dcs_s
            ws.cell(r, 15, f"=IF(N{r}=0,NA(),100*M{r}/N{r})")  # pct
            ws.cell(r, 16, f"=OR(E{r}>{pulse_threshold},F{r}>{pulse_threshold})")  # is_pulse

        ws.column_dimensions["A"].width = 22
        for col in range(2, len(headers) + 1):
            ws.column_dimensions[get_column_letter(col)].width = 18
        n = len(times) + 1
        return n

    last_rows: Dict[str, int] = {}
    for i, (side, title) in enumerate(zip(sides, titles)):
        last_rows[title[:31]] = _write_side_sheet(title, i)

    ws_sum = wb.create_sheet("Summary", 1)
    ws_sum["A1"] = "Non-pulse aggregation (Excel SUMIF on is_pulse = FALSE)"
    ws_sum["A1"].font = Font(bold=True)
    ws_sum["A2"] = (
        "Balance check: sum_dcs_actual should nearly equal sum_dcs_attributed; "
        "sum_dcs_residual is the finite-step leftover. "
        "pct_abs_uptake_vs_abs_s = 100 * |sum dcs_uptake| / |sum dcs_s|."
    )
    sum_headers = [
        "side",
        "sum_dcs_actual_mol_m3",
        "sum_dcs_uptake_mol_m3",
        "sum_dcs_s_mol_m3",
        "sum_dcs_attributed_mol_m3",
        "sum_dcs_residual_mol_m3",
        "pct_abs_uptake_vs_abs_s",
        "sum_dMS_mol",
    ]
    for col, name in enumerate(sum_headers, start=1):
        cell = ws_sum.cell(4, col, name)
        cell.font = Font(bold=True)

    for row_i, title in enumerate(titles, start=5):
        sheet = title[:31]
        n = last_rows[sheet]
        # is_pulse is column P; attributed terms H/I/J/K/L; dMS_mol is G.
        ws_sum.cell(row_i, 1, sheet)
        ws_sum.cell(row_i, 2, f"=SUMIF('{sheet}'!P3:P{n},FALSE,'{sheet}'!H3:H{n})")
        ws_sum.cell(row_i, 3, f"=SUMIF('{sheet}'!P3:P{n},FALSE,'{sheet}'!I3:I{n})")
        ws_sum.cell(row_i, 4, f"=SUMIF('{sheet}'!P3:P{n},FALSE,'{sheet}'!J3:J{n})")
        ws_sum.cell(row_i, 5, f"=SUMIF('{sheet}'!P3:P{n},FALSE,'{sheet}'!K3:K{n})")
        ws_sum.cell(row_i, 6, f"=SUMIF('{sheet}'!P3:P{n},FALSE,'{sheet}'!L3:L{n})")
        ws_sum.cell(row_i, 7, f"=IF(D{row_i}=0,NA(),100*ABS(C{row_i})/ABS(D{row_i}))")
        ws_sum.cell(row_i, 8, f"=SUMIF('{sheet}'!P3:P{n},FALSE,'{sheet}'!G3:G{n})")

    for col, width in zip("ABCDEFGH", [28, 22, 22, 20, 24, 22, 24, 14]):
        ws_sum.column_dimensions[col].width = width

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "salt_concentration_attribution.xlsx"
    wb.save(out_path)
    return out_path


# Config ``plots`` key -> plotting function.
PLOT_FUNCTIONS: Dict[str, Callable[..., Path]] = {
    "sap_flux_vs_uptake": plot_sap_flux_vs_uptake,
    "daily_uptake_vs_sap_flux": plot_daily_uptake_vs_sap_flux,
    "relative_water_content": plot_relative_water_content,
    "soil_salt": plot_soil_salt,
    "salt_concentrations": plot_salt_concentrations,
    "matric_potential": plot_matric_potential,
}

# Config ``plots`` key -> non-figure export (Excel / tables).
EXPORT_FUNCTIONS: Dict[str, Callable[..., Path]] = {
    "salt_concentration_attribution": export_salt_concentration_attribution,
}


def main(argv: Optional[List[str]] = None) -> List[Path]:
    """Make every enabled figure (and configured exports) for a batch.

    Returns:
        Paths of the saved figures and export files.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="Path to the lab simulation config.yaml")
    parser.add_argument("--batch", default="latest", help="Batch id under results/, or 'latest'")
    parser.add_argument("--only", nargs="+", help="Subset of plot/export keys to make")
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
        if key in EXPORT_FUNCTIONS:
            path = EXPORT_FUNCTIONS[key](df, cfg, plot_cfg, batch_dir / "figures")
        elif key in PLOT_FUNCTIONS:
            path = PLOT_FUNCTIONS[key](df, cfg, plot_cfg, common, batch_dir / "figures")
        else:
            raise KeyError(f"No plot/export function registered for '{key}'")
        print(f"  Saved {path.relative_to(batch_dir)}")
        saved.append(path)
    return saved


if __name__ == "__main__":
    main()
