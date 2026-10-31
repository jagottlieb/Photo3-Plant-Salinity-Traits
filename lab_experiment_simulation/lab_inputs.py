"""
Lab data inputs for the lab experiment simulation.

Turns the processed lab data into model inputs: effective compartment depth,
initial relative soil moisture and salinity, irrigation and salt pulses, and
sap flux on a ground-area basis for comparison with root water uptake.

Units: water fluxes in um/s (model ``qs``), salt mass in mol/m2 of ground
(model ``SaltySoilMultiple.MS``), concentrations in mM (= mol/m3).

Effective rooting depth
-----------------------
The soil model keeps salt as an area density ``MS`` (mol/m2) and updates
concentration as ``cs = MS / (s * N * ZR)``. Lab pots have a measured soil
volume ``V`` and ground area ``A``, so I set an effective depth

    ZR_eff = V / A

per compartment. Then ``cs = MS / (s * N * ZR_eff)`` is the same as
``cs = M_total / (V * N * s)`` with the pot moles, without changing the
ground-area bookkeeping used by uptake and plant salt pools. Sensitivity
runs keep a uniform species rooting depth (0.2 m for Pecan) instead.
"""

from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
CM2_TO_M2 = 1e-4
ML_TO_M3 = 1e-6
L_TO_M3 = 1e-3


def _read_csv(lab: Dict[str, Any], key: str, time_col: str) -> pd.DataFrame:
    """Read a lab data CSV from ``lab[key]`` (repo-relative), indexed by its datetime column."""
    return pd.read_csv(REPO_ROOT / lab[key], parse_dates=[time_col]).set_index(time_col).sort_index()


def initial_cutoff(lab: Dict[str, Any], date: pd.Timestamp) -> pd.Timestamp | None:
    """First water or salt pulse on ``date``, when initials should stop before it.

    With ``lab["initial_before_pulses"]`` set, the start-date mean uses only samples
    from before this time so a same-day addition is not already inside the initial state.
    """
    if not lab.get("initial_before_pulses"):
        return None
    df = _read_csv(lab, "irrigation_file", "timestamp")
    day = df[df.index.normalize() == pd.Timestamp(date).normalize()]
    prefix = f"tree{lab['tree']}"
    amount_cols = [f"{prefix}_{side}_{kind}" for side in lab["sides"] for kind in ("vol_liters", "moles")]
    active = day[day[amount_cols].fillna(0).sum(axis=1) > 0]
    if active.empty:
        return None
    return pd.Timestamp(active.index[0])


def _day_mean(series: pd.Series, date: pd.Timestamp, label: str, before: pd.Timestamp | None = None) -> float:
    """Mean of ``series`` over one calendar day; error if that window is empty.

    ``before`` drops samples at or after that time, used when a pulse falls on the start date.
    """
    day = series[series.index.normalize() == pd.Timestamp(date).normalize()].dropna()
    if before is not None:
        day = day[day.index < before]
    if day.empty:
        when = f" before {before}" if before is not None else ""
        raise ValueError(f"No {label} data on {pd.Timestamp(date).date()}{when}")
    return float(day.mean())


def container_volumes_mL(lab: Dict[str, Any]) -> List[float]:
    """Soil volume (mL) per compartment from ``lab["container_volume_mL"]``.

    Accepts a per-side dict (e.g. ``{left: 2600, right: 3900}``) or a scalar
    applied to every side.
    """
    vol = lab["container_volume_mL"]
    if isinstance(vol, dict):
        return [float(vol[side]) for side in lab["sides"]]
    return [float(vol)] * len(lab["sides"])


def effective_zr_m(lab: Dict[str, Any]) -> List[float]:
    """Effective compartment depths ZR_eff = V / A_ground (m), one per side.

    Matches pot salt mass balance to the model's ``MS`` (mol/m2) form:
    ``cs = MS / (s * N * ZR_eff)`` equals ``M_total / (V * N * s)``.
    """
    area_m2 = lab["ground_area_cm2"] * CM2_TO_M2
    return [vol_mL * ML_TO_M3 / area_m2 for vol_mL in container_volumes_mL(lab)]


def compartment_depth(lab: Dict[str, Any]) -> float:
    """Scalar ZR_eff (m) when ``container_volume_mL`` is a single value.

    Prefer ``effective_zr_m`` for lab runs (supports per-side volumes).
    """
    depths = effective_zr_m(lab)
    if len(set(round(z, 12) for z in depths)) != 1:
        raise TypeError(
            "container_volume_mL differs by side; use effective_zr_m() for ZR_eff"
        )
    return depths[0]


def _vwc_depths_for_side(lab: Dict[str, Any], side: str) -> List[str]:
    """Probe depths used for initial s on one side.

    ``lab["vwc_depths"]`` may be a shared list (e.g. ``[shallow]``), a per-side
    dict (e.g. ``{left: deep, right: shallow}`` or ``{left: [deep], ...}``), or
    omitted (default both shallow and deep).
    """
    depths = lab.get("vwc_depths", ["shallow", "deep"])
    if isinstance(depths, dict):
        side_depths = depths[side]
        if isinstance(side_depths, str):
            return [side_depths]
        return list(side_depths)
    return list(depths)


def initial_s(lab: Dict[str, Any], date: pd.Timestamp, porosity: float) -> List[float]:
    """Relative soil moisture per compartment from the start-date mean VWC.

    Each side's VWC is the mean of the depths from ``_vwc_depths_for_side``,
    averaged over ``date`` and divided by ``porosity``. With
    ``initial_before_pulses``, the mean stops before the first pulse that day.
    """
    df = _read_csv(lab, "vwc_file", "datetime")
    before = initial_cutoff(lab, date)
    s = []
    for side in lab["sides"]:
        depths = _vwc_depths_for_side(lab, side)
        cols = [f"tree{lab['tree']}_{side}_{depth}" for depth in depths]
        vwc = _day_mean(df[cols].mean(axis=1), date, f"{side} VWC", before)
        s.append(vwc / porosity)
    return s


def initial_cs(lab: Dict[str, Any], date: pd.Timestamp) -> List[float]:
    """Soil salinity (mM) per compartment: start-date mean of each side's conductivity probe.

    With ``initial_before_pulses``, the mean stops before the first pulse that day.
    """
    df = _read_csv(lab, "conductivity_file", "datetime")
    before = initial_cutoff(lab, date)
    return [
        _day_mean(df[lab["conductivity_cols"][side]], date, f"{side} conductivity", before)
        for side in lab["sides"]
    ]


def pulse_schedule(lab: Dict[str, Any], start: pd.Timestamp, end: pd.Timestamp, porosity: float,
                   step_times: pd.DatetimeIndex) -> Dict[int, Dict[str, np.ndarray]]:
    """Irrigation and salt pulses between ``start`` and ``end``, keyed by model step.

    Each pulse is applied at the start of the first model step beginning at or
    after the irrigation timestamp. Water is converted to relative moisture as
    ``V_irr / (V_container * porosity)``; salt to mol per m2 of ground.

    Args:
        lab: The config ``lab`` section.
        start, end: Window in which pulses are applied.
        porosity: Soil porosity N (-).
        step_times: Start time of each model step.

    Returns:
        ``{step: {"ds": [...], "dms": [...], "time": Timestamp}}`` with one value per compartment.
    """
    df = _read_csv(lab, "irrigation_file", "timestamp")
    df = df[(df.index >= start) & (df.index < end)]
    prefix = f"tree{lab['tree']}"
    vol = df[[f"{prefix}_{side}_vol_liters" for side in lab["sides"]]].to_numpy(dtype=float)
    mol = df[[f"{prefix}_{side}_moles" for side in lab["sides"]]].to_numpy(dtype=float)
    pore_volume = np.array(container_volumes_mL(lab), dtype=float) * ML_TO_M3 * porosity
    ground_area = lab["ground_area_cm2"] * CM2_TO_M2

    pulses = {}
    for t, v, m in zip(df.index, vol, mol):
        if not (v.any() or m.any()):
            continue
        step = int(np.searchsorted(step_times, t))
        pulses[step] = {"ds": v * L_TO_M3 / pore_volume, "dms": m / ground_area, "time": t}
    return pulses


def apply_pulse(soil: Any, ds: np.ndarray, dms: np.ndarray) -> None:
    """Add water (relative moisture) and salt (mol/m2) to each soil compartment.

    Moisture is capped at saturation; all added salt stays in the soil.
    """
    soil.s = np.minimum(soil.s + ds, 1.0)
    soil.MS = soil.MS + dms
    soil.cs = soil.MS / (soil.s * soil.N * soil.ZR)


def sap_flux_ground_basis(lab: Dict[str, Any]) -> pd.DataFrame:
    """Sap flux per side converted to a ground-area basis (um/s, same units as ``qs``).

    Sap velocity (mm/s) is multiplied by ``root_area_cm2`` / ground area.
    """
    df = _read_csv(lab, "sap_flux_file", "datetime")
    out = pd.DataFrame(index=df.index)
    areas = lab.get("root_area_cm2")
    if not areas:
        raise KeyError("lab.root_area_cm2 is required (cm2 per side)")
    for side in lab["sides"]:
        sapwood_cm2 = areas[side]
        out[side] = df[f"tree{lab['tree']}_{side}_vs"] * sapwood_cm2 / lab["ground_area_cm2"] * 1000
    return out
