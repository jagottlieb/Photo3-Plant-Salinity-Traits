#!/usr/bin/env python
"""Compare Berger.leak() to Hyprop unsaturated conductivity measurements.

Two figures share the same curves and points:
    berger_k_vs_theta.png  x = volumetric water content
    berger_k_vs_s.png      x = relative soil moisture s = theta / N

Curves:
    Berger.leak (Mualem–van Genuchten, parameters already in soil.py)
    Clapp–Hornberger fit, K = KS * s^(2B+3), with N fixed at Berger.N
    Empirical fit, log10 K = a + b * theta

Both fits are ordinary least squares on log10 K. The evaluation sheet
``Evaluation-Conductivity K(pF)`` stores tension (kPa) and log10 K (cm/d).
Those rows are the same points as ``Evaluation-Conductivity K(Θ)``, which
carries the paired water content (identical log10 K, same order).

``Berger.leak`` returns K in µm/s (cm/d times 0.11574). The curve is converted
back to cm/d so it sits on the Hyprop axis.

Usage (from the repo root):
    python lab_experiment_simulation/diagnostics/conductivity_curve/plot_berger_k_theta.py
"""

import sys
from pathlib import Path
from typing import NamedTuple, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

DIAG_DIR = Path(__file__).resolve().parent
REPO_ROOT = DIAG_DIR.parents[2]
sys.path.insert(0, str(REPO_ROOT))

from soil import Berger

HYPROP_XLSX = (
    REPO_ROOT
    / "data"
    / "soil_data"
    / "hyprop_data_and_reten_curves"
    / "Berger_Sample_points_removed.xlsx"
)
KPF_SHEET = "Evaluation-Conductivity K(pF)"
# soil.Berger.leak multiplies K [cm/d] by this factor and returns µm/s.
CM_D_TO_UM_S = 0.11574
OUT_THETA = DIAG_DIR / "berger_k_vs_theta.png"
OUT_S = DIAG_DIR / "berger_k_vs_s.png"


def set_plot_style(font_size: float = 11) -> None:
    """Apply the project plot style (seaborn whitegrid, serif font stack)."""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Caladea", "STIX Two Text", "Constantia", "Georgia", "Cambria"],
        "font.size": font_size,
        "mathtext.fontset": "stix",
    })


def load_conductivity_points(path: Path = HYPROP_XLSX) -> pd.DataFrame:
    """Hyprop evaluation points: tension and log10 K from the K(pF) sheet, theta from K(Θ).

    Args:
        path: Hyprop workbook (Berger sample, points removed).

    Returns:
        One row per evaluation point, with theta (m3/m3) and K (cm/d).
    """
    xl = pd.ExcelFile(path)
    k_pf = pd.read_excel(xl, KPF_SHEET)
    k_theta_name = next(
        name for name in xl.sheet_names
        if name.startswith("Evaluation-Conductivity K") and "pF" not in name
    )
    k_theta = pd.read_excel(xl, k_theta_name)
    log10_k = k_pf["log 10 K [cm/d]"].to_numpy(dtype=float)
    log10_k_theta = k_theta["log 10 K [cm/d]"].to_numpy(dtype=float)
    if log10_k.shape != log10_k_theta.shape or not np.allclose(log10_k, log10_k_theta):
        raise ValueError(f"{KPF_SHEET} and {k_theta_name} log10 K columns do not match")
    theta = k_theta["Water Content [Vol%]"].to_numpy(dtype=float) / 100.0
    return pd.DataFrame({
        "tension_kpa": k_pf["[kPa]"].to_numpy(dtype=float),
        "theta": theta,
        "log10_k_cm_d": log10_k,
        "k_cm_d": 10.0 ** log10_k,
    })


class ClappHornbergerFit(NamedTuple):
    """K = ks * s^(2*b + 3), fitted in log10 K with porosity fixed."""

    ks: float
    b: float
    rmse_log10: float


class LogThetaFit(NamedTuple):
    """log10 K = intercept + slope * theta, fitted by least squares."""

    intercept: float
    slope: float
    rmse_log10: float


def rmse_log10(observed: np.ndarray, predicted_log10: np.ndarray) -> float:
    """Root-mean-square error of log10 K."""
    return float(np.sqrt(np.mean((observed - predicted_log10) ** 2)))


def fit_clapp_hornberger(s: np.ndarray, log10_k: np.ndarray) -> ClappHornbergerFit:
    """Fit K = KS * s^(2B+3) by linear regression of log10 K on log10 s.

    Porosity is not a free parameter: s is already theta / N.
    """
    exponent, log10_ks = np.polyfit(np.log10(s), log10_k, 1)
    ks = float(10.0 ** log10_ks)
    b = float((exponent - 3.0) / 2.0)
    predicted = np.log10(clapp_hornberger_k(s, ks, b))
    return ClappHornbergerFit(ks, b, rmse_log10(log10_k, predicted))


def clapp_hornberger_k(s: np.ndarray, ks: float, b: float) -> np.ndarray:
    """Clapp–Hornberger conductivity (cm/d) at relative soil moisture s."""
    return ks * np.asarray(s, dtype=float) ** (2.0 * b + 3.0)


def fit_log_theta(theta: np.ndarray, log10_k: np.ndarray) -> LogThetaFit:
    """Fit log10 K = a + b*theta by linear regression."""
    slope, intercept = np.polyfit(theta, log10_k, 1)
    predicted = intercept + slope * theta
    return LogThetaFit(float(intercept), float(slope), rmse_log10(log10_k, predicted))


def log_theta_k(theta: np.ndarray, fit: LogThetaFit) -> np.ndarray:
    """Conductivity (cm/d) from the log-linear water-content fit."""
    return 10.0 ** (fit.intercept + fit.slope * np.asarray(theta, dtype=float))


def model_k_cm_d(theta: np.ndarray, texture: Berger) -> np.ndarray:
    """Unsaturated K (cm/d) from Berger.leak at volumetric water content theta.

    Relative saturation passed to leak is s = theta / N, which leak converts
    back to theta before the Mualem-van Genuchten calculation.
    """
    s = np.asarray(theta, dtype=float) / texture.N
    return np.asarray(texture.leak(s), dtype=float) / CM_D_TO_UM_S


def plot_k(
    x_line: np.ndarray,
    curves: Sequence[Tuple[np.ndarray, str, str]],
    x_points: np.ndarray,
    k_points: np.ndarray,
    xlabel: str,
    title: str,
    out_path: Path,
) -> plt.Figure:
    """Plot conductivity curves and the Hyprop points against one moisture x-axis."""
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for k_line, label, color in curves:
        ax.plot(x_line, k_line, color=color, lw=2, label=label, zorder=2)
    ax.scatter(
        x_points,
        k_points,
        s=36,
        color="#d95f0e",
        edgecolor="white",
        linewidth=0.6,
        label="Hyprop evaluation points",
        zorder=3,
    )
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Hydraulic conductivity, $K$ (cm d$^{-1}$)")
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper left", fontsize="small")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")
    return fig


def plot_conductivity(points: pd.DataFrame, texture: Berger) -> None:
    """Fit both alternatives and save the water-content and relative-moisture figures."""
    pad = 0.008
    theta = points["theta"].to_numpy()
    log10_k = points["log10_k_cm_d"].to_numpy()
    theta_line = np.linspace(theta.min() - pad, theta.max() + pad, 400)
    s = theta / texture.N
    s_line = theta_line / texture.N
    ch = fit_clapp_hornberger(s, log10_k)
    log_fit = fit_log_theta(theta, log10_k)
    mualem_log10 = np.log10(model_k_cm_d(theta, texture))

    curves_theta = [
        (model_k_cm_d(theta_line, texture), "Berger.leak (Mualem–van Genuchten)", "#08519c"),
        (clapp_hornberger_k(s_line, ch.ks, ch.b), fr"Clapp–Hornberger ($K_S$={ch.ks:.0f}, $B$={ch.b:.2f})", "#31a354"),
        (log_theta_k(theta_line, log_fit), rf"$\log_{{10}} K = {log_fit.intercept:.2f} + {log_fit.slope:.2f}\,\theta$", "#6a51a3"),
    ]
    plot_k(
        theta_line,
        curves_theta,
        theta,
        points["k_cm_d"].to_numpy(),
        r"Volumetric water content, $\theta$ (m$^3$ m$^{-3}$)",
        "Berger soil: hydraulic conductivity vs water content",
        OUT_THETA,
    )
    plot_k(
        s_line,
        curves_theta,
        s,
        points["k_cm_d"].to_numpy(),
        rf"Relative soil moisture, $s = \theta / N$ ($N$ = {texture.N:g})",
        "Berger soil: hydraulic conductivity vs relative soil moisture",
        OUT_S,
    )

    print(f"Clapp-Hornberger, N = {texture.N:g} fixed: K = {ch.ks:.1f} * s^({2 * ch.b + 3:.2f})")
    print(f"  KS = {ch.ks:.1f} cm/d, B = {ch.b:.3f}, RMSE(log10 K) = {ch.rmse_log10:.3f}")
    print(f"log10 K = {log_fit.intercept:.3f} + {log_fit.slope:.3f} * theta")
    print(f"  RMSE(log10 K) = {log_fit.rmse_log10:.3f}")
    print(f"Mualem-van Genuchten RMSE(log10 K) = {rmse_log10(log10_k, mualem_log10):.3f}")
    print(f"Hyprop relative soil moisture range: s = {s.min():.3f} to {s.max():.3f}")


if __name__ == "__main__":
    set_plot_style()
    texture = Berger()
    plot_conductivity(load_conductivity_points(), texture)
