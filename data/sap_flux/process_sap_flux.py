"""Load, concatenate and tidy hourly sap velocity data for Tree 1.

Conversions:
    "Left_Vs_out_mm_s" / "Right_Vs_out_mm_s" [mm/s] -> tree1_left_vs / tree1_right_vs [mm/s]  (rename only)
"""
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DATA_DIR = Path(__file__).parent
RAW_GLOB = "sap_velocity_hourly_mm_s_*.csv"

TREE = 1
SIDES = ["left", "right"]


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
    """Read and concatenate all hourly sap velocity .csv files, indexed by datetime.

    Files are ordered by their first timestamp (filenames are inconsistent), and
    duplicate timestamps at export boundaries keep the row from the later file.
    """
    frames = []
    for f in data_dir.glob(RAW_GLOB):
        df = pd.read_csv(f, parse_dates=["datetime"])
        df["source_file"] = f.name
        frames.append(df)
    frames.sort(key=lambda d: d["datetime"].min())
    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates("datetime", keep="last")
    return df.set_index("datetime").sort_index()


def convert(df: pd.DataFrame, tree: int = TREE) -> pd.DataFrame:
    """Rename "Side_Vs_out_mm_s" columns to treeN_side_vs (velocity stays in mm/s)."""
    rename = {f"{s.title()}_Vs_out_mm_s": f"tree{tree}_{s}_vs" for s in SIDES}
    return df.rename(columns=rename)


def vs_cols(tree: int = TREE) -> list[str]:
    """Sap velocity column names for one tree."""
    return [f"tree{tree}_{s}_vs" for s in SIDES]


def plot(df: pd.DataFrame, tree: int = TREE, out_dir: Path = DATA_DIR) -> plt.Figure:
    """Plot left and right sap velocity timeseries for one tree and save as PNG.

    Missing hours are left as gaps rather than interpolated across.
    """
    file_breaks = df.groupby("source_file").apply(lambda g: g.index.min()).sort_values().iloc[1:]
    df = df.asfreq("h")
    fig, ax = plt.subplots(figsize=(12, 4))

    for col in vs_cols(tree):
        ax.plot(df.index, df[col], label=col)
    for t in file_breaks:
        ax.axvline(t, color="grey", ls=":", lw=1)

    ax.set_ylabel("Sap velocity (mm/s)")
    ax.legend(loc="best", fontsize="small", frameon=False)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.set_xlabel("Date (2026)")
    fig.suptitle(f"Sap velocity (Tree {tree})")
    fig.tight_layout()
    path = out_dir / f"sap_velocity_tree{tree}_timeseries.png"
    fig.savefig(path, dpi=150)
    print(f"Saved {path}")
    return fig


if __name__ == "__main__":
    set_plot_style()
    df = convert(load_raw())
    df[vs_cols()].round(8).to_csv(DATA_DIR / f"sap_velocity_tree{TREE}_processed.csv")
    plot(df)
    plt.show()
