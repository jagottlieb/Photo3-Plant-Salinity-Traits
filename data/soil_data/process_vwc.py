"""Load, concatenate and tidy Berger hourly volumetric water content data.

Deep probes (and trees 2-4 shallow) come from the recalibrated CSV exports:
    "Tree N - Side Depth" [m3/m3] -> treeN_side_depth [m3/m3]  (rename only)

Tree-1 shallow probes are rebuilt from Board ID 9 TEROS raw voltage (mV) in
``matric_potential_conductivity/`` (soil2 = left, soil3 = right), using the
Probe-2 cubic calibration:

    VWC = 2.138593e-9*mV^3 - 1.486886e-5*mV^2 + 3.461576e-2*mV - 2.669084e1

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

from process_matric_conductivity import _read_excel

DATA_DIR = Path(__file__).parent / "vwc_data"
MATRIC_DIR = Path(__file__).parent / "matric_potential_conductivity"
ARCHIVE_DIR = DATA_DIR / "archive"

TREES = [1, 2, 3, 4]
SIDES = ["left", "right"]
DEPTHS = ["shallow", "deep"]
PLOT_TREE = 1

# Board ID 9 TEROS raw-voltage columns -> tree-1 shallow sides (same mapping as conductivity).
SHALLOW_MV_COLS = {"left": "soil2_vwc", "right": "soil3_vwc"}

# Probe-2 cubic: mV -> volumetric water content (m3/m3).
SHALLOW_VWC_CAL = (2.138593e-9, -1.486886e-5, 3.461576e-2, -2.669084e1)

# Tree 1 pot geometry (matches lab_experiment_simulation configs).
TREE1_VOLUME_ML = {"left": 2600.0, "right": 3900.0}
TREE1_GROUND_AREA_CM2 = 210.0
ML_TO_M3 = 1e-6
CM2_TO_M2 = 1e-4

PROCESSED_ARTIFACTS = (
    "vwc_processed.csv",
    "vwc_flux.csv",
    "vwc_tree1_timeseries.png",
)


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
    """Read and concatenate all hourly recalibrated .csv files, indexed by datetime.

    Exports overlap and the final hour of an earlier export can be a partial
    average, so duplicate timestamps keep the row from the latest file.
    """
    frames = []
    for f in sorted(data_dir.glob("Berger*.csv")):
        if "recalibrated" not in f.name:
            continue
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


def mv_to_vwc(mv: pd.Series) -> pd.Series:
    """Convert TEROS raw voltage (mV) to VWC with the Probe-2 cubic calibration."""
    a, b, c, d = SHALLOW_VWC_CAL
    return a * mv**3 + b * mv**2 + c * mv + d


def load_shallow_from_mv(matric_dir: Path = MATRIC_DIR) -> pd.DataFrame:
    """Tree-1 left/right shallow VWC from Board ID 9 TEROS raw voltage timeseries.

    ``soil2_vwc`` / ``soil3_vwc`` are raw mV (same series used for pore-water EC).
    Duplicate timestamps keep the row from the latest export file.
    """
    frames = []
    for f in sorted(matric_dir.glob("Berger Hourly Data Collection_Board_ID 9_*.xlsx")):
        df = _read_excel(f)
        df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
        keep = ["datetime", *SHALLOW_MV_COLS.values()]
        df = df[keep].copy()
        df["source_file_mv"] = f.name
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No Board ID 9 workbooks in {matric_dir}")
    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates("datetime", keep="last")
    df = df.set_index("datetime").sort_index()
    out = pd.DataFrame(index=df.index)
    out["source_file_mv"] = df["source_file_mv"]
    for side, col in SHALLOW_MV_COLS.items():
        out[f"tree1_{side}_shallow"] = mv_to_vwc(df[col])
        out[f"tree1_{side}_shallow_mv"] = df[col]
    return out


def apply_shallow_mv_calibration(df: pd.DataFrame, shallow: pd.DataFrame) -> pd.DataFrame:
    """Replace tree-1 shallow columns with mV-calibrated values; union the indexes."""
    out = df.reindex(df.index.union(shallow.index)).sort_index()
    for side in SIDES:
        col = f"tree1_{side}_shallow"
        out[col] = shallow[col].reindex(out.index)
        out[f"{col}_mv"] = shallow[f"{col}_mv"].reindex(out.index)
    if "source_file_mv" in shallow.columns:
        out["source_file_mv"] = shallow["source_file_mv"].reindex(out.index)
    return out


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
    src = df["source_file"] if "source_file" in df.columns else df.get("source_file_mv")
    file_breaks = (
        df.assign(_src=src).dropna(subset=["_src"])
        .groupby("_src").apply(lambda g: g.index.min(), include_groups=False)
        .sort_values().iloc[1:]
        if src is not None else pd.Series(dtype="datetime64[ns]")
    )
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


def plot_shallow_calibration_compare(
    old: pd.DataFrame,
    new: pd.DataFrame,
    out_dir: Path = DATA_DIR,
) -> plt.Figure:
    """Compare archived vs mV-calibrated tree-1 left/right shallow VWC timeseries."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    for ax, side in zip(axes, SIDES):
        col = f"tree1_{side}_shallow"
        old_s = old[col].asfreq("h") if col in old.columns else None
        new_s = new[col].asfreq("h")
        if old_s is not None:
            ax.plot(old_s.index, old_s, color="0.55", lw=1.2, label=f"{side} shallow (old recalibrated)")
        ax.plot(new_s.index, new_s, color="C0" if side == "left" else "C1",
                lw=1.2, label=f"{side} shallow (mV cubic)")
        ax.set_ylabel(f"{side.title()} shallow VWC\n(m$^3$/m$^3$)")
        ax.legend(loc="best", fontsize="small", frameon=False)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    axes[-1].set_xlabel("Date (2026)")
    fig.suptitle("Tree 1 shallow VWC: old recalibrated vs Board-9 mV cubic")
    fig.tight_layout()
    path = out_dir / "vwc_tree1_shallow_calibration_compare.png"
    fig.savefig(path, dpi=150)
    print(f"Saved {path}")
    return fig


def archive_processed(data_dir: Path = DATA_DIR, archive_dir: Path = ARCHIVE_DIR) -> list[Path]:
    """Move current processed VWC artifacts into archive/ before regenerating."""
    archive_dir.mkdir(parents=True, exist_ok=True)
    moved: list[Path] = []
    for name in PROCESSED_ARTIFACTS:
        src = data_dir / name
        if not src.exists():
            continue
        dest = archive_dir / name
        if dest.exists():
            dest.unlink()
        src.replace(dest)
        moved.append(dest)
        print(f"Archived {src.name} -> {dest}")
    return moved


def build() -> pd.DataFrame:
    """Load recalibrated deep/other probes and overlay tree-1 shallow from Board-9 mV."""
    df = convert(load_raw())
    shallow = load_shallow_from_mv()
    return apply_shallow_mv_calibration(df, shallow)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")

    set_plot_style()
    old_path = ARCHIVE_DIR / "vwc_processed.csv"
    live_path = DATA_DIR / "vwc_processed.csv"
    if live_path.exists() and not old_path.exists():
        # Snapshot current products first so the diagnostic can load the prior shallow series.
        archive_processed()

    old = (
        pd.read_csv(old_path, parse_dates=["datetime"]).set_index("datetime").sort_index()
        if old_path.exists() else None
    )

    df = build()
    df[vwc_cols()].round(6).to_csv(DATA_DIR / "vwc_processed.csv")
    print(f"Saved {DATA_DIR / 'vwc_processed.csv'}", flush=True)
    flux = vwc_flux(df, tree=1)
    flux.round(12).to_csv(DATA_DIR / "vwc_flux.csv")
    print(f"Saved {DATA_DIR / 'vwc_flux.csv'}", flush=True)
    plot(df)
    if old is not None:
        plot_shallow_calibration_compare(old, df)
