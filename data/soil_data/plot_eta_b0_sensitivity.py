"""Compare Hilhorst pore-water EC at several dielectric offsets.

eta_b0 is the bulk permittivity at zero moisture, subtracted in the denominator:

    cond_pw = eta_p * bulk_EC / (eta_b - eta_b0)

Usage (from the repo root):
    python data/soil_data/plot_eta_b0_sensitivity.py
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_matric_conductivity as proc

OFFSETS = (4.1, 3.0, 2.0, 1.7)
COLORS = {
    4.1: "#b2182b",
    3.0: "#ef8a62",
    2.0: "#67a9cf",
    1.7: "#2166ac",
}
OUT = proc.DATA_DIR / "porewater_eta_b0_sensitivity.png"
# soil2 at 4.1 reaches hundreds of dS/m where eta_b approaches the offset.
YLIM = (0, 4.5)


def main() -> None:
    proc.set_plot_style()
    df = proc.load_raw()
    file_breaks = df.groupby("source_file").apply(lambda g: g.index.min(), include_groups=False).sort_values().iloc[1:]
    hourly = df.asfreq("h")

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, sharey=True)
    for ax, col in zip(axes, proc.CONDUCT_COLS):
        num = col.removeprefix("soil").split("_", 1)[0]
        bulk = hourly[col] / 1000
        raw = hourly[f"soil{num}_vwc"]
        temp = hourly[f"soil{num}_temp"]
        for offset in OFFSETS:
            series = proc.pore_water_ec(bulk, raw, temp, eta_b0=offset)
            label = rf"$\eta_{{b0}} = {offset:g}$"
            if offset == proc.ETA_B0:
                label += " (current)"
            ax.plot(hourly.index, series, color=COLORS[offset], lw=1.2, label=label)
        for t in file_breaks:
            ax.axvline(t, color="grey", ls=":", lw=1)
        ax.set_ylabel("Pore-water EC (dS/m)")
        ax.set_ylim(*YLIM)
        ax.set_title(f"soil{num}")
        ax.legend(loc="upper right", fontsize="small", frameon=False, ncol=2)

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    axes[-1].set_xlabel("Date (2026)")
    fig.suptitle(
        r"Hilhorst pore-water EC at several $\eta_{b0}$"
        "\n"
        r"soil2 at $\eta_{b0}=4.1$ leaves this scale (peak 368 dS/m) where $\eta_b$ nears 4.1"
    )
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
