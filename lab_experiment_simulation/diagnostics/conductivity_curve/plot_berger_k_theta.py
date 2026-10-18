#!/usr/bin/env python
"""Compare Berger.leak() to Hyprop unsaturated conductivity measurements.

Two figures share the same curve and points:
    berger_k_vs_theta.png  x = volumetric water content
    berger_k_vs_s.png      x = relative soil moisture s = theta / N

The evaluation sheet ``Evaluation-Conductivity K(pF)`` stores tension (kPa) and
log10 K (cm/d). Those rows are the same points as ``Evaluation-Conductivity K(Θ)``,
which carries the paired water content (identical log10 K, same order).

``Berger.leak`` returns K in µm/s (cm/d times 0.11574). The curve is converted
back to cm/d so it sits on the Hyprop axis.

Usage (from the repo root):
    python lab_experiment_simulation/diagnostics/conductivity_curve/plot_berger_k_theta.py
"""

import sys
from pathlib import Path

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


def model_k_cm_d(theta: np.ndarray, texture: Berger) -> np.ndarray:
    """Unsaturated K (cm/d) from Berger.leak at volumetric water content theta.

    Relative saturation passed to leak is s = theta / N, which leak converts
    back to theta before the Mualem-van Genuchten calculation.
    """
    s = np.asarray(theta, dtype=float) / texture.N
    return np.asarray(texture.leak(s), dtype=float) / CM_D_TO_UM_S


def plot_k(
    x_line: np.ndarray,
    k_line: np.ndarray,
    x_points: np.ndarray,
    k_points: np.ndarray,
    xlabel: str,
    title: str,
    out_path: Path,
) -> plt.Figure:
    """Plot Berger.leak and the Hyprop points against one moisture x-axis."""
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(
        x_line,
        k_line,
        color="#08519c",
        lw=2,
        label="Berger.leak (Mualem–van Genuchten)",
        zorder=2,
    )
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
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    return fig


def plot_conductivity(points: pd.DataFrame, texture: Berger) -> None:
    """Save the volumetric-water-content figure and the relative-moisture figure."""
    pad = 0.008
    theta = points["theta"].to_numpy()
    theta_line = np.linspace(theta.min() - pad, theta.max() + pad, 400)
    s = theta / texture.N
    s_line = theta_line / texture.N
    k_line = model_k_cm_d(theta_line, texture)
    k_points = points["k_cm_d"].to_numpy()

    plot_k(
        theta_line,
        k_line,
        theta,
        k_points,
        r"Volumetric water content, $\theta$ (m$^3$ m$^{-3}$)",
        "Berger soil: hydraulic conductivity vs water content",
        OUT_THETA,
    )
    plot_k(
        s_line,
        k_line,
        s,
        k_points,
        rf"Relative soil moisture, $s = \theta / N$ ($N$ = {texture.N:g})",
        "Berger soil: hydraulic conductivity vs relative soil moisture",
        OUT_S,
    )

    comparison = points.copy()
    comparison["s"] = s
    comparison["log10_k_model"] = np.log10(model_k_cm_d(theta, texture))
    comparison["model_over_data"] = 10.0 ** comparison["log10_k_model"] / points["k_cm_d"].to_numpy()
    printable = comparison[["theta", "s", "log10_k_cm_d", "log10_k_model", "model_over_data"]]
    print(printable.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"Hyprop relative soil moisture range: s = {s.min():.3f} to {s.max():.3f}")


if __name__ == "__main__":
    set_plot_style()
    texture = Berger()
    plot_conductivity(load_conductivity_points(), texture)
