"""Compare no-salt + no-evap against salt-treatment + no-evap.

Figures (written into the diagnostic batch ``figures/``):
  qs_hr_compare.png             per side: qs, cumulative-HR bars, then one psi_s panel
                                per run with matric and osmotic fills
  qs_hr_compare_total_psi.png   same, with a single panel of total psi_s for both runs

From the repo root:
  python lab_experiment_simulation/salt_treatment/diagnostics/no_salt_additions_no_evap/plot_compare.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parents[4]
DIAG_DIR = Path(__file__).resolve().parent
# Matching no-evap salt run (overlap-day initials); resolved to latest completed batch.
BASELINE_CONFIG_DIR = (
    REPO_ROOT / "lab_experiment_simulation" / "salt_treatment" / "diagnostics" / "no_evap"
)
sys.path.insert(0, str(REPO_ROOT / "sensitivity_analysis"))
sys.path.append(str(REPO_ROOT / "lab_experiment_simulation"))
sys.path.append(str(REPO_ROOT))

from post_processing import set_plot_style  # sensitivity_analysis/post_processing.py
import soil

# Lab pulse annotation / legend helpers (same boxes as sap_flux_vs_uptake).
_lab_pp_spec = importlib.util.spec_from_file_location(
    "lab_post_processing",
    REPO_ROOT / "lab_experiment_simulation" / "post_processing.py",
)
_lab_pp = importlib.util.module_from_spec(_lab_pp_spec)
assert _lab_pp_spec.loader is not None
_lab_pp_spec.loader.exec_module(_lab_pp)

BASELINE_COLOR = "#08519c"
NOSALT_COLOR = "#6baed6"  # light blue, distinct from official dark blue
MATRIC_COLOR = "#6baed6"
OSMOTIC_COLOR = "#969696"
WITH_SALT_LABEL = "Model: With salt addition"
NO_SALT_LABEL = "Model: No salt addition"

# qs_hr_compare layout. Panel heights are in units of PANEL_UNIT_IN; the
# potential panels are 30% shorter than the qs panel.
QS_H, HR_H, PSI_H = 2.2, 1.0, 2.2 * 0.70
PANEL_UNIT_IN = 1.04   # inches per height unit
PANEL_GAP_IN = 0.14    # between panels of one side
SIDE_GAP_IN = 0.45     # between the left-side and right-side blocks
# Cumulative-HR axis, shared by both sides (mm).
HR_LIMIT_MM, HR_TICK_MM = 0.2, 0.05


def _savefig(fig: plt.Figure, path: Path, dpi: int) -> Path:
    """Save a figure, using a temp file when OneDrive/IDE has the target locked."""
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    except OSError:
        tmp = path.with_suffix(path.suffix + ".tmp.png")
        fig.savefig(tmp, dpi=dpi, bbox_inches="tight")
        tmp.replace(path)
    print(f"Saved {path}")
    return path


def _load_cfg(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _resolve_latest_batch(config_dir: Path) -> Path:
    results_dir = config_dir / "results"
    batches = sorted(p for p in results_dir.iterdir() if any((p / "timeseries").glob("*.csv")))
    if not batches:
        raise FileNotFoundError(f"No completed batches in {results_dir}")
    return batches[-1]


def _load_timeseries(batch_dir: Path, name: str) -> pd.DataFrame:
    df = pd.read_csv(batch_dir / "timeseries" / f"{name}.csv", parse_dates=["datetime"])
    return df[~df["is_burn_in"]].set_index("datetime").sort_index()


def _compartment(df: pd.DataFrame, variable: str, i: int) -> np.ndarray:
    return np.array([json.loads(cell)[i] for cell in df[variable]], dtype=float)


def _dt_seconds(df: pd.DataFrame) -> np.ndarray:
    """Per-row timestep length (s). First step uses the median of later steps."""
    if "time_days" in df.columns:
        t = df["time_days"].to_numpy(dtype=float)
        dt_s = np.diff(t, prepend=np.nan) * 86400.0
    else:
        dt_s = df.index.to_series().diff().dt.total_seconds().to_numpy(dtype=float)
    finite = dt_s[np.isfinite(dt_s) & (dt_s > 0)]
    fill = float(np.median(finite)) if len(finite) else 3600.0
    if not np.isfinite(dt_s[0]) or dt_s[0] <= 0:
        dt_s[0] = fill
    return dt_s


def _cumulative_hr_mm(df: pd.DataFrame, i: int) -> pd.Series:
    """Cumulative hydraulic redistribution into compartment ``i`` (mm).

    HR rate is ``-min(qs, 0)`` (µm/s). Integrate with ``um/s * s / 1000 -> mm``.
    """
    qs = _compartment(df, "qs", i)
    hr = np.where(qs < 0.0, -qs, 0.0)
    return pd.Series(np.cumsum(hr * _dt_seconds(df) / 1000.0), index=df.index)


def _daily_cumulative_hr(hr_cum: pd.Series) -> pd.Series:
    """End-of-day cumulative HR (mm), indexed by calendar day at noon."""
    if hr_cum.empty:
        return pd.Series(dtype=float)
    daily = hr_cum.groupby(hr_cum.index.normalize()).last()
    daily.index = daily.index + pd.Timedelta(hours=12)
    return daily


def _annotate_side_pulses(ax: plt.Axes, irrigation: pd.DataFrame, lab: Dict[str, Any],
                          side: str, placement: str, pulse_y: float) -> Tuple[bool, bool]:
    """Draw irrigation/salt pulse markers; return (has_volume, has_salt)."""
    prefix = f"tree{lab['tree']}_{side}"
    has_volume = has_salt = False
    for t, row in irrigation[irrigation[f"{prefix}_vol_liters"] > 0].iterrows():
        note = f"{row[f'{prefix}_vol_liters']} L"
        has_volume = True
        # Label experiment inputs from the irrigation sheet (official salt moles).
        if row[f"{prefix}_moles"] > 0:
            note += f", {row[f'{prefix}_moles']} mol"
            has_salt = True
        _lab_pp._annotate_pulse(ax, t, note, placement=placement, y=pulse_y)
    return has_volume, has_salt


def _pulse_vlines(ax: plt.Axes, irrigation: pd.DataFrame, lab: Dict[str, Any],
                  side: str) -> None:
    """Mark irrigation times on an axis without labels (for the HR panel)."""
    prefix = f"tree{lab['tree']}_{side}"
    for t in irrigation.index[irrigation[f"{prefix}_vol_liters"] > 0]:
        ax.axvline(t, color="0.5", ls=":", lw=0.8, zorder=0)


def _soil_texture(cfg: Dict[str, Any]) -> Any:
    """Soil texture object with the config's retention-curve formulation."""
    sim_cfg = cfg["simulation"]
    texture = getattr(soil, sim_cfg["soil_texture"])()
    formulation = sim_cfg.get("soil_formulation")
    if formulation is not None and hasattr(texture, "FORMULATION"):
        texture.FORMULATION = str(formulation).strip().lower().replace("-", "_")
    return texture


def _psi_components(df: pd.DataFrame, i: int, texture: Any
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Matric, osmotic, and total soil water potential (MPa) for compartment ``i``."""
    s = _compartment(df, "s", i)
    matric = np.asarray(texture.psi_s(s), dtype=float)
    total = _compartment(df, "psi_s", i)
    osmotic = total - matric
    return matric, osmotic, total


def _plot_psi_stack(ax: plt.Axes, df: pd.DataFrame, i: int, texture: Any,
                    linewidth: float, legend: bool) -> None:
    """Stacked matric + osmotic fill under total ``psi_s``."""
    matric, _osmotic, total = _psi_components(df, i, texture)
    # zorder 1.5 keeps the fills above the daily lines from _apply_date_axis (zorder 1).
    ax.fill_between(df.index, 0.0, matric, color=MATRIC_COLOR, alpha=0.85,
                    label="Matric", linewidth=0, zorder=1.5)
    ax.fill_between(df.index, matric, total, color=OSMOTIC_COLOR, alpha=0.85,
                    label="Osmotic", linewidth=0, zorder=1.5)
    ax.plot(df.index, total, color="0.15", lw=linewidth, label=r"Total $\psi_s$", zorder=2)
    ax.axhline(0, color="0.7", lw=0.8, zorder=0)
    if legend:
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize="small", frameon=False)


def _plot_psi_total(ax: plt.Axes, baseline: pd.DataFrame, nosalt: pd.DataFrame, i: int,
                    linewidth: float) -> None:
    """Total ``psi_s`` for both runs, in the same colors as the qs panel."""
    ax.plot(baseline.index, _compartment(baseline, "psi_s", i), color=BASELINE_COLOR,
            lw=linewidth, label=WITH_SALT_LABEL, zorder=2)
    ax.plot(nosalt.index, _compartment(nosalt, "psi_s", i), color=NOSALT_COLOR,
            lw=linewidth, label=NO_SALT_LABEL, zorder=2)
    ax.axhline(0, color="0.7", lw=0.8, zorder=0)


def _side_grids(n_sides: int, height_ratios: List[float], width: float
                ) -> Tuple[plt.Figure, List[Any]]:
    """Figure plus one stacked sub-grid per side, with panel sizes and gaps fixed in inches.

    The figure height follows the panel count, so adding or dropping a panel
    changes the figure size rather than the size of the other panels.
    """
    n = len(height_ratios)
    panels_in = PANEL_UNIT_IN * sum(height_ratios)
    block_in = panels_in + (n - 1) * PANEL_GAP_IN
    usable_in = n_sides * block_in + (n_sides - 1) * SIDE_GAP_IN
    top, bottom = plt.rcParams["figure.subplot.top"], plt.rcParams["figure.subplot.bottom"]
    fig = plt.figure(figsize=(width, usable_in / (top - bottom)))
    # GridSpec hspace is a fraction of the average cell height.
    outer = GridSpec(n_sides, 1, figure=fig, hspace=SIDE_GAP_IN / block_in)
    inner_hspace = PANEL_GAP_IN / (panels_in / n)
    grids = [
        outer[k].subgridspec(n, 1, height_ratios=height_ratios, hspace=inner_hspace)
        for k in range(n_sides)
    ]
    return fig, grids


def plot_qs_hr_compare(baseline: pd.DataFrame, nosalt: pd.DataFrame, cfg: Dict[str, Any],
                       out_dir: Path, psi_mode: str = "components") -> Path:
    """Per side: qs, downward cumulative-HR bars, then soil water potential.

    Args:
        psi_mode: ``components`` gives one panel per run with matric and osmotic
            fills under total psi_s (``qs_hr_compare.png``). ``total`` gives a
            single panel per side with both runs' total psi_s
            (``qs_hr_compare_total_psi.png``).
    """
    if psi_mode not in ("components", "total"):
        raise ValueError(f"Unknown psi_mode '{psi_mode}'")
    components = psi_mode == "components"
    lab = cfg["lab"]
    plot_cfg = cfg.get("plots", {}).get("sap_flux_vs_uptake", {})
    common = cfg.get("plots", {}).get("common", {})
    irrigation = _lab_pp._load_irrigation(lab, baseline.index[0], baseline.index[-1])
    texture = _soil_texture(cfg)
    linewidth = float(common.get("linewidth", 1.2))
    width, _height = common.get("figsize", [9, 3.0])
    placement = str(plot_cfg.get("pulse_placement", "in_pane"))
    pulse_y = float(plot_cfg.get("pulse_y", 0.92))
    bar_width = 0.35
    # With-salt left of noon, no-salt right.
    offsets = {WITH_SALT_LABEL: -0.2, NO_SALT_LABEL: 0.2}
    n_psi = 2 if components else 1
    height_ratios = [QS_H, HR_H] + [PSI_H] * n_psi

    fig, grids = _side_grids(len(lab["sides"]), height_ratios, width)
    titles = ["Left side", "Right side"]
    psi_run_labels = (WITH_SALT_LABEL, NO_SALT_LABEL)
    psi_runs = (baseline, nosalt)
    bottom_ax = None

    # Each row uses the same y extents on the left and right sides, taken from
    # both runs and both sides.
    n_sides = len(lab["sides"])
    qs_all = [_compartment(df, "qs", k) for df in psi_runs for k in range(n_sides)]
    qs_lo = float(min(0.0, min(v.min() for v in qs_all)))
    qs_hi = float(max(v.max() for v in qs_all))
    qs_span = max(qs_hi - qs_lo, 1e-6)
    qs_ylim = (qs_lo - 0.05 * qs_span, qs_hi + 0.28 * qs_span)  # headroom for pulse labels
    psi_all = [_compartment(df, "psi_s", k) for df in psi_runs for k in range(n_sides)]
    psi_lo = float(min(v.min() for v in psi_all))
    psi_hi = float(max(0.0, max(v.max() for v in psi_all)))
    psi_span = max(psi_hi - psi_lo, 1e-6)
    psi_ylim = (psi_lo - 0.05 * psi_span, psi_hi + 0.08 * psi_span)

    for i, side in enumerate(lab["sides"]):
        inner = grids[i]
        ax_qs = fig.add_subplot(inner[0], sharex=bottom_ax)
        ax_hr = fig.add_subplot(inner[1], sharex=ax_qs)
        ax_psis = [fig.add_subplot(inner[2 + j], sharex=ax_qs) for j in range(n_psi)]
        bottom_ax = ax_psis[-1]

        qs_base = _compartment(baseline, "qs", i)
        qs_alt = _compartment(nosalt, "qs", i)
        ax_qs.plot(baseline.index, qs_base, color=BASELINE_COLOR,
                   lw=linewidth, label=WITH_SALT_LABEL)
        ax_qs.plot(nosalt.index, qs_alt, color=NOSALT_COLOR,
                   lw=linewidth, label=NO_SALT_LABEL)
        ax_qs.axhline(0, color="0.7", lw=0.8)
        ax_qs.set_ylim(qs_ylim)

        has_volume, has_salt = _annotate_side_pulses(
            ax_qs, irrigation, lab, side, placement, pulse_y,
        )
        ax_qs.set_title(titles[i], fontsize="medium")
        ax_qs.set_ylabel(r"Root uptake $q_s$ ($\mu$m s$^{-1}$)")
        _lab_pp._legend(
            ax_qs, common,
            **_lab_pp._pulse_legend_kwargs(_lab_pp._pulse_legend_text(has_volume, has_salt)),
        )
        plt.setp(ax_qs.get_xticklabels(), visible=False)

        # Daily end-of-day cumulative HR as bars growing downward from zero.
        series = (
            (WITH_SALT_LABEL, BASELINE_COLOR, _daily_cumulative_hr(_cumulative_hr_mm(baseline, i))),
            (NO_SALT_LABEL, NOSALT_COLOR, _daily_cumulative_hr(_cumulative_hr_mm(nosalt, i))),
        )
        for label, color, daily in series:
            if daily.empty:
                continue
            x = mdates.date2num(daily.index.to_list()) + offsets[label]
            ax_hr.bar(
                x, -daily.to_numpy(), width=bar_width, align="center",
                color=color, edgecolor=color, linewidth=0.6, alpha=0.9,
                label=label, zorder=2,
            )
        _pulse_vlines(ax_hr, irrigation, lab, side)
        ax_hr.axhline(0, color="0.7", lw=0.8, zorder=1)
        # Keep zero at the top; magnitude increases downward. Tick labels stay positive mm.
        ax_hr.set_ylim(-HR_LIMIT_MM, 0.0)
        ax_hr.set_yticks(-np.arange(0.0, HR_LIMIT_MM + 0.5 * HR_TICK_MM, HR_TICK_MM))
        ax_hr.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda y, _pos: f"{abs(y):g}")
        )
        ax_hr.set_ylabel("Cumulative HR (mm)")
        plt.setp(ax_hr.get_xticklabels(), visible=False)

        if components:
            for j, (ax_psi, df, run_label) in enumerate(zip(ax_psis, psi_runs, psi_run_labels)):
                _plot_psi_stack(
                    ax_psi, df, i, texture, linewidth,
                    legend=(i == 0 and j == 0),
                )
                ax_psi.set_ylabel(f"{run_label}\n" + r"$\psi_s$ (MPa)")
        else:
            _plot_psi_total(ax_psis[0], baseline, nosalt, i, linewidth)
            ax_psis[0].set_ylabel(r"Soil $\psi_s$ (MPa)")
        for j, ax_psi in enumerate(ax_psis):
            _pulse_vlines(ax_psi, irrigation, lab, side)
            ax_psi.set_ylim(psi_ylim)
            if not (i == len(lab["sides"]) - 1 and j == len(ax_psis) - 1):
                plt.setp(ax_psi.get_xticklabels(), visible=False)

    _lab_pp._apply_date_axis(bottom_ax, common)
    bottom_ax.set_xlabel("Date (2026)")

    stem = "qs_hr_compare" if components else "qs_hr_compare_total_psi"
    path = _savefig(fig, out_dir / f"{stem}.png", int(common.get("dpi", 300)))
    plt.close(fig)
    return path


def main() -> None:
    set_plot_style()
    cfg = _load_cfg(DIAG_DIR / "config.yaml")
    diag_batch = _resolve_latest_batch(DIAG_DIR)
    baseline_batch = _resolve_latest_batch(BASELINE_CONFIG_DIR)
    name = cfg["simulation"]["name"]
    baseline = _load_timeseries(baseline_batch, name)
    nosalt = _load_timeseries(diag_batch, name)
    out_dir = diag_batch / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Baseline (salt, no evap): {baseline_batch}")
    print(f"Diagnostic (no salt, no evap): {diag_batch}")
    plot_qs_hr_compare(baseline, nosalt, cfg, out_dir, psi_mode="components")
    plot_qs_hr_compare(baseline, nosalt, cfg, out_dir, psi_mode="total")


if __name__ == "__main__":
    main()
