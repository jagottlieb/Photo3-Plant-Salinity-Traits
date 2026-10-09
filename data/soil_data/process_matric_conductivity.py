"""Load, concatenate and convert Berger hourly matric potential / conductivity data.

Conversions:
    soilN_matric  [kPa]   -> soilN_matric_MPa  [MPa]   (/ 1000)
    soilN_conduct [uS/cm] -> soilN_conduct_dSm [dS/m]  (/ 1000)
    soilN_conduct_dSm     -> soilN_conduct_mM  [mM]    (* 10)
"""
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DATA_DIR = Path(__file__).parent / "matric_potential_conductivity"

MATRIC_COLS = ["soil0_matric", "soil1_matric"]
CONDUCT_COLS = ["soil2_conduct", "soil3_conduct"]


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
    """Read and concatenate all hourly .xlsx files, indexed by datetime."""
    frames = []
    for f in sorted(data_dir.glob("*.xlsx")):
        df = pd.read_excel(f)
        df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
        df["source_file"] = f.name
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df = df.sort_values("datetime").drop_duplicates("datetime", keep="first")
    return df.set_index("datetime")


def convert(df: pd.DataFrame) -> pd.DataFrame:
    """Add matric potential in MPa and conductivity in dS/m and mM."""
    df = df.copy()
    for col in MATRIC_COLS:
        df[f"{col}_MPa"] = df[col] / 1000
    for col in CONDUCT_COLS:
        df[f"{col}_dSm"] = df[col] / 1000
        df[f"{col}_mM"] = df[f"{col}_dSm"] * 10
    return df


def plot(df: pd.DataFrame, out_dir: Path = DATA_DIR) -> plt.Figure:
    """Plot matric potential and conductivity timeseries and save as PNG."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)

    for col in MATRIC_COLS:
        axes[0].plot(df.index, df[f"{col}_MPa"], label=col)
    axes[0].set_ylabel("Matric potential (MPa)")

    for col in CONDUCT_COLS:
        axes[1].plot(df.index, df[f"{col}_dSm"], label=col)
    axes[1].set_ylabel("Conductivity (dS/m)")

    for col in CONDUCT_COLS:
        axes[2].plot(df.index, df[f"{col}_mM"], label=col)
    axes[2].set_ylabel("Conductivity (mM)")

    file_breaks = df.groupby("source_file").apply(lambda g: g.index.min()).iloc[1:]
    for ax in axes:
        for t in file_breaks:
            ax.axvline(t, color="grey", ls=":", lw=1)
        ax.legend(loc="best", fontsize="small", frameon=False)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    axes[-1].set_xlabel("Date (2026)")
    fig.suptitle("Berger soil: matric potential and conductivity (Board ID 9)")
    fig.tight_layout()
    path = out_dir / "matric_conductivity_timeseries.png"
    fig.savefig(path, dpi=150)
    print(f"Saved {path}")
    return fig


if __name__ == "__main__":
    set_plot_style()
    df = convert(load_raw())
    df.to_csv(DATA_DIR / "matric_conductivity_processed.csv")
    plot(df)
    plt.show()
