"""Load, concatenate and tidy Berger hourly recalibrated volumetric water content data.

Conversions:
    "Tree N - Side Depth" [m3/m3] -> treeN_side_depth [m3/m3]  (rename only)

Tree-1 VWC change is also turned into a ground-area flux (um/s), depletion positive:
    vwc_difference          = theta(t+1) - theta(t)
    vol_flux                = -vwc_difference * V_side / dt     [m3/s]
    vol_flux_scaled_to_root = vol_flux / A_ground * 1e6        [um/s]
"""
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DATA_DIR = Path(__file__).parent / "vwc_data"

TREES = [1, 2, 3, 4]
SIDES = ["left", "right"]
DEPTHS = ["shallow", "deep"]
PLOT_TREE = 1

# Tree 1 pot geometry (matches lab_experiment_simulation configs).
TREE1_VOLUME_ML = {"left": 2600.0, "right": 3900.0}
TREE1_GROUND_AREA_CM2 = 210.0
ML_TO_M3 = 1e-6
CM2_TO_M2 = 1e-4


def set_plot_style(font_size: float = 10) -> None:
    """Apply the project plot style (seaborn whitegrid, serif font stack)."""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Caladea", "STIX Two Text", "Constantia", "Georgia", "Cambria"],
        "font.size": font_size,
        "mathtext.fontset": "stix",
    })


def load_raw(data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Read and concatenate all hourly .csv files, indexed by datetime.

    Exports overlap and the final hour of an earlier export can be a partial
    average, so duplicate timestamps keep the row from the latest file.
    """
    frames = []
    for f in sorted(data_dir.glob("Berger*.csv")):
        df = pd.read_csv(f, parse_dates=["Time"])
        df["source_file"] = f.name
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df = df.rename(columns={"Time": "datetime"})
    df = df.drop_duplicates("datetime", keep="last")
    return df.set_index("datetime").sort_index()


def convert(df: pd.DataFrame) -> pd.DataFrame:
    """Rename "Tree N - Side Depth" columns to treeN_side_depth (VWC stays in m3/m3)."""
    df = df.copy()
    rename = {
        f"Tree {t} - {s.title()} {d.title()}": f"tree{t}_{s}_{d}"
        for t in TREES for s in SIDES for d in DEPTHS
    }
    return df.rename(columns=rename)


def vwc_cols(tree: int | None = None) -> list[str]:
    """VWC column names, optionally restricted to a single tree."""
    trees = TREES if tree is None else [tree]
    return [f"tree{t}_{s}_{d}" for t in trees for d in DEPTHS for s in SIDES]


def vwc_flux(df: pd.DataFrame, tree: int = 1) -> pd.DataFrame:
    """Forward VWC difference to ground-area flux (um/s); depletion is positive.

    Uses tree-1 pot volumes and ground area. Columns per probe:
    ``{col}_vwc_difference``, ``{col}_vol_flux`` (m3/s),
    ``{col}_vol_flux_scaled_to_root`` (um/s).
    """
    if tree != 1:
        raise ValueError("vwc_flux geometry is only defined for tree 1")
    area_m2 = TREE1_GROUND_AREA_CM2 * CM2_TO_M2
    out = pd.DataFrame(index=df.index)
    dt_s = df.index.to_series().diff(-1).dt.total_seconds().abs()
    for side in SIDES:
        vol_m3 = TREE1_VOLUME_ML[side] * ML_TO_M3
        for depth in DEPTHS:
            col = f"tree{tree}_{side}_{depth}"
            dtheta = df[col].shift(-1) - df[col]
            vol = -dtheta * vol_m3 / dt_s
            out[f"{col}_vwc_difference"] = dtheta
            out[f"{col}_vol_flux"] = vol
            out[f"{col}_vol_flux_scaled_to_root"] = vol / area_m2 * 1e6
    return out


def plot(df: pd.DataFrame, tree: int = PLOT_TREE, out_dir: Path = DATA_DIR) -> plt.Figure:
    """Plot shallow and deep VWC timeseries for one tree and save as PNG.

    Missing hours are left as gaps rather than interpolated across.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    file_breaks = df.groupby("source_file").apply(lambda g: g.index.min()).sort_values().iloc[1:]
    df = df.asfreq("h")
    fig, axes = plt.subplots(len(DEPTHS), 1, figsize=(12, 6), sharex=True, sharey=True)

    for ax, depth in zip(axes, DEPTHS):
        for side in SIDES:
            col = f"tree{tree}_{side}_{depth}"
            ax.plot(df.index, df[col], label=col)
        ax.set_ylabel(f"{depth.title()} VWC (m$^3$/m$^3$)")

    for ax in axes:
        for t in file_breaks:
            ax.axvline(t, color="grey", ls=":", lw=1)
        ax.legend(loc="best", fontsize="small", frameon=False)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    axes[-1].set_xlabel("Date (2026)")
    fig.suptitle(f"Berger soil: volumetric water content (Tree {tree})")
    fig.tight_layout()
    path = out_dir / f"vwc_tree{tree}_timeseries.png"
    fig.savefig(path, dpi=150)
    print(f"Saved {path}")
    return fig


if __name__ == "__main__":
    set_plot_style()
    df = convert(load_raw())
    df[vwc_cols()].round(6).to_csv(DATA_DIR / "vwc_processed.csv")
    flux = vwc_flux(df, tree=1)
    flux.round(12).to_csv(DATA_DIR / "vwc_flux.csv")
    print(f"Saved {DATA_DIR / 'vwc_flux.csv'}")
    plot(df)
    plt.show()
