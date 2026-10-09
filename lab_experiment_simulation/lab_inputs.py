"""
Lab data inputs for the lab experiment simulation.

Turns the processed lab data into model inputs: compartment depth, initial
relative soil moisture and salinity, irrigation and salt pulses, and sap flux
on a ground-area basis for comparison with root water uptake.

Units: water fluxes in um/s (model ``qs``), salt mass in mol/m2 of ground
(model ``SaltySoilMultiple.MS``), concentrations in mM (= mol/m3).
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


def _day_mean(series: pd.Series, date: pd.Timestamp, label: str) -> float:
    """Mean of ``series`` over one calendar day; error if there is no data that day."""
    day = series[series.index.normalize() == date.normalize()].dropna()
    if day.empty:
        raise ValueError(f"No {label} data on {date.date()}")
    return float(day.mean())


def compartment_depth(lab: Dict[str, Any]) -> float:
    """Soil depth (m) of one compartment: container volume / ground area."""
    return lab["container_volume_mL"] * ML_TO_M3 / (lab["ground_area_cm2"] * CM2_TO_M2)


def initial_s(lab: Dict[str, Any], date: pd.Timestamp, porosity: float) -> List[float]:
    """Relative soil moisture per compartment from the start-date mean VWC.

    Each side's VWC is the mean of its shallow and deep probes (or the one probe
    with data), averaged over ``date`` and divided by ``porosity``.
    """
    df = _read_csv(lab, "vwc_file", "datetime")
    s = []
    for side in lab["sides"]:
        cols = [f"tree{lab['tree']}_{side}_{depth}" for depth in ("shallow", "deep")]
        vwc = _day_mean(df[cols].mean(axis=1), date, f"{side} VWC")
        s.append(vwc / porosity)
    return s


def initial_cs(lab: Dict[str, Any], date: pd.Timestamp) -> List[float]:
    """Soil salinity (mM) per compartment: start-date mean of each side's conductivity probe."""
    df = _read_csv(lab, "conductivity_file", "datetime")
    return [_day_mean(df[lab["conductivity_cols"][side]], date, f"{side} conductivity") for side in lab["sides"]]


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
    pore_volume = lab["container_volume_mL"] * ML_TO_M3 * porosity
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

    Sap velocity (mm/s) is multiplied by root sapwood area / ground area, with
    the whole root cross-section (pi d^2 / 4) treated as sapwood.
    """
    df = _read_csv(lab, "sap_flux_file", "datetime")
    out = pd.DataFrame(index=df.index)
    for side in lab["sides"]:
        sapwood_cm2 = np.pi * lab["root_diameter_cm"][side] ** 2 / 4
        out[side] = df[f"tree{lab['tree']}_{side}_vs"] * sapwood_cm2 / lab["ground_area_cm2"] * 1000
    return out
