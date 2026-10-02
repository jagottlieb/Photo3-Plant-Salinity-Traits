"""Schematic of the dynamic filtration efficiency (FE) adjustment used in hydraulics.py.

FE is computed by calling the model's own uptake() methods, so the curves stay in
sync with any revision of the dynamic FE rule. The activation curve replaces the
model's pre-threshold baseline with a lower, illustrative value (E0_ACTIVATION).
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from hydraulics import HalophyteStemLeafStorageMultiComp, HalophyteStemStorageMultiComp  # noqa: E402

E0 = 0.5  # baseline filtration efficiency (-)
E0_ACTIVATION = 0.3  # pre-threshold FE for the activation curve (-)
C_MAX = 80.0  # concentration at which FE reaches its upper limit (mM)
E_UPPER = 0.9999  # cap on FE in uptake()
C_RANGE = (0.0, 100.0)

TICK_LABELS = {
    "c_thresh": "Threshold\nconcentration",
    "c_max": "Maximum FE\nconcentration",
    "E0": "Baseline FE",
    "E0_activation": "Pre-activation FE",
    "E_upper": "Maximum FE",
}

# (class, line style, legend label, pre-threshold FE override or None)
CURVES = [
    (HalophyteStemStorageMultiComp, "k-", "Continuous increase in efficiency", None),
    (
        HalophyteStemLeafStorageMultiComp,
        "k--",
        "Activation of filtration then continuous increase in efficiency",
        E0_ACTIVATION,
    ),
]

SAVE_DIR = Path(__file__).resolve().parent / "Figures"


def dynamic_fe(cls: type, c_arr: np.ndarray, e0: float, c_max: float) -> np.ndarray:
    """Back out the dynamic FE (-) from cls.uptake() for each stem concentration in c_arr."""
    obj = cls.__new__(cls)
    obj.E = e0
    obj.dynamic_E = True
    obj.c_stem_max = c_max
    obj.Salt_Uptake = True
    fe = np.empty_like(c_arr, dtype=float)
    for i, c in enumerate(c_arr):
        obj.c_stem = c
        u = obj.uptake([1.0], [1.0])
        u0 = obj.uptake([1.0], [1.0], E=0.0)
        fe[i] = 1.0 - u / u0
    return fe


def set_plot_style() -> None:
    """Match the photosynthesis reduction figure style."""
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "font.size": 11,
        "axes.titlesize": 11,
        "axes.labelsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
    })


def main() -> None:
    """Plot FE versus salt concentration for both dynamic FE formulations."""
    set_plot_style()
    c_arr = np.linspace(*C_RANGE, 500)
    c_thresh = E0 * C_MAX

    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    for cls, style, label, pre_thresh_fe in CURVES:
        fe = dynamic_fe(cls, c_arr, E0, C_MAX)
        if pre_thresh_fe is not None:
            fe[c_arr < c_thresh] = pre_thresh_fe
        ax.plot(c_arr, fe, style, label=label)

    guide = dict(color="gray", linestyle=":", linewidth=1, zorder=0)
    for x in (c_thresh, C_MAX):
        ax.axvline(x, **guide)
    for y in (E0_ACTIVATION, E0, E_UPPER):
        ax.axhline(y, **guide)

    ax.set_xlim(*C_RANGE)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([c_thresh, C_MAX], labels=[TICK_LABELS["c_thresh"], TICK_LABELS["c_max"]])
    ax.set_yticks(
        [E0_ACTIVATION, E0, E_UPPER],
        labels=[TICK_LABELS["E0_activation"], TICK_LABELS["E0"], TICK_LABELS["E_upper"]],
    )
    ax.minorticks_off()
    ax.spines[["top", "right"]].set_visible(False)

    ax.set_xlabel("Salt Concentration (mMol)")
    ax.set_ylabel("Filtration Efficiency (-)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False)

    plt.tight_layout()
    SAVE_DIR.mkdir(exist_ok=True)
    fig.savefig(SAVE_DIR / "dynamic_filtration_efficiency.png", dpi=300, bbox_inches="tight")
    fig.savefig(SAVE_DIR / "dynamic_filtration_efficiency.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
