"""Load, concatenate and convert Berger hourly matric potential / conductivity data.

Conversions:
    soilN_matric  [kPa]   -> soilN_matric_MPa  [MPa]   (/ 1000)
    soilN_conduct [uS/cm] -> soilN_conduct_dSm [dS/m]  (/ 1000)
    soilN_conduct_dSm     -> soilN_conduct_mM  [mM]    (* 10)
    soilN_conduct_dSm, soilN_vwc (TEROS raw count), soilN_temp [C]
                          -> cond_pw_N         [dS/m]  (Hilhorst pore-water EC)
    cond_pw_N             -> cond_pw_N_mM      [mM]    (* 10)

Pore-water EC uses the TEROS 11/12 mineral-soil calibration. soilN_vwc is the
raw count (RAW), not a water fraction. eta_b0 = 4.1.

    eta_p = 80.3 - 0.37 * (T_soil - 20)
    eta_b = (2.887e-9*RAW^3 - 2.080e-5*RAW^2 + 5.276e-2*RAW - 43.39)^2
    cond_pw = eta_p * bulk_EC / (eta_b - eta_b0)
"""
import ctypes
import io
from ctypes import wintypes
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DATA_DIR = Path(__file__).parent / "matric_potential_conductivity"

MATRIC_COLS = ["soil0_matric", "soil1_matric"]
CONDUCT_COLS = ["soil2_conduct", "soil3_conduct"]
ETA_B0 = 4.1  # dielectric offset in the Hilhorst pore-water model


def pore_water_column(conduct_col: str) -> str:
    """cond_pw_N for a soilN_conduct column."""
    num = conduct_col.removeprefix("soil").split("_", 1)[0]
    return f"cond_pw_{num}"


def set_plot_style(font_size: float = 10) -> None:
    """Apply the project plot style (seaborn whitegrid, serif font stack)."""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Caladea", "STIX Two Text", "Constantia", "Georgia", "Cambria"],
        "font.size": font_size,
        "mathtext.fontset": "stix",
    })


def _read_excel(path: Path) -> pd.DataFrame:
    """Read an hourly workbook, including one another process still has open."""
    try:
        return pd.read_excel(path)
    except PermissionError:
        return pd.read_excel(io.BytesIO(_read_shared(path)))


def _read_shared(path: Path) -> bytes:
    """Read a file Windows has locked against a normal open (Excel or OneDrive)."""
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.CreateFileW.argtypes = [
        wintypes.LPCWSTR, wintypes.DWORD, wintypes.DWORD, wintypes.LPVOID,
        wintypes.DWORD, wintypes.DWORD, wintypes.HANDLE,
    ]
    kernel.CreateFileW.restype = wintypes.HANDLE
    kernel.GetFileSizeEx.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.LARGE_INTEGER)]
    kernel.GetFileSizeEx.restype = wintypes.BOOL
    kernel.ReadFile.argtypes = [
        wintypes.HANDLE, wintypes.LPVOID, wintypes.DWORD,
        ctypes.POINTER(wintypes.DWORD), wintypes.LPVOID,
    ]
    kernel.ReadFile.restype = wintypes.BOOL
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel.CloseHandle.restype = wintypes.BOOL
    handle = kernel.CreateFileW(
        str(path.resolve()), 0x80000000, 0x1 | 0x2 | 0x4, None, 3, 0x80, None,
    )
    if handle == wintypes.HANDLE(-1).value:
        raise PermissionError(f"Cannot read {path.name} (winerror {ctypes.get_last_error()})")
    try:
        size = wintypes.LARGE_INTEGER()
        if not kernel.GetFileSizeEx(handle, ctypes.byref(size)):
            raise PermissionError(f"Cannot size {path.name}")
        buf = (ctypes.c_char * size.value)()
        read = wintypes.DWORD()
        if not kernel.ReadFile(handle, buf, size.value, ctypes.byref(read), None):
            raise PermissionError(f"Cannot read {path.name}")
        return bytes(buf[: read.value])
    finally:
        kernel.CloseHandle(handle)


def load_raw(data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Read and concatenate all hourly .xlsx files, indexed by datetime.

    Exports overlap and the final hour of an earlier export can be a partial
    average, so duplicate timestamps keep the row from the latest file.
    """
    frames = []
    for f in sorted(data_dir.glob("*.xlsx")):
        df = _read_excel(f)
        df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
        df["source_file"] = f.name
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates("datetime", keep="last")
    return df.set_index("datetime").sort_index()


def dielectric_from_raw(raw: pd.Series) -> pd.Series:
    """Bulk dielectric permittivity from TEROS raw counts (mineral-soil calibration)."""
    inner = 2.887e-9 * raw**3 - 2.080e-5 * raw**2 + 5.276e-2 * raw - 43.39
    return inner**2


def pore_water_ec(bulk_dsm: pd.Series, raw: pd.Series, temp_c: pd.Series) -> pd.Series:
    """Pore-water EC (dS/m) from bulk EC, TEROS raw counts, and soil temperature (C).

    Rows with eta_b <= eta_b0 are left blank; the Hilhorst denominator is not positive there.
    """
    eta_p = 80.3 - 0.37 * (temp_c - 20.0)
    denom = dielectric_from_raw(raw) - ETA_B0
    return (eta_p * bulk_dsm / denom).where(denom > 0)


def convert(df: pd.DataFrame) -> pd.DataFrame:
    """Add matric potential in MPa, bulk conductivity in dS/m and mM, and pore-water EC in dS/m and mM."""
    df = df.copy()
    for col in MATRIC_COLS:
        df[f"{col}_MPa"] = df[col] / 1000
    for col in CONDUCT_COLS:
        num = col.removeprefix("soil").split("_", 1)[0]
        df[f"{col}_dSm"] = df[col] / 1000
        df[f"{col}_mM"] = df[f"{col}_dSm"] * 10
        pw = pore_water_column(col)
        df[pw] = pore_water_ec(df[f"{col}_dSm"], df[f"soil{num}_vwc"], df[f"soil{num}_temp"])
        df[f"{pw}_mM"] = df[pw] * 10
    return df


def plot(df: pd.DataFrame, out_dir: Path = DATA_DIR) -> plt.Figure:
    """Plot matric potential and conductivity timeseries and save as PNG.

    Missing hours are left as gaps rather than interpolated across.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    file_breaks = df.groupby("source_file").apply(lambda g: g.index.min()).sort_values().iloc[1:]
    df = df.asfreq("h")
    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

    for col in MATRIC_COLS:
        axes[0].plot(df.index, df[f"{col}_MPa"], label=col)
    axes[0].set_ylabel("Matric potential (MPa)")

    for col in CONDUCT_COLS:
        axes[1].plot(df.index, df[f"{col}_dSm"], label=col)
    axes[1].set_ylabel("Bulk EC (dS/m)")

    for col in CONDUCT_COLS:
        axes[2].plot(df.index, df[pore_water_column(col)], label=pore_water_column(col))
    axes[2].set_ylabel("Pore-water EC (dS/m)")
    axes[2].set_ylim(0, 10)  # hide the dry-end Hilhorst spikes while we inspect the series

    for col in CONDUCT_COLS:
        axes[3].plot(df.index, df[f"{col}_mM"], label=col)
    axes[3].set_ylabel("Bulk EC (mM)")

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
    pw_cols = [name for c in CONDUCT_COLS for name in (pore_water_column(c), f"{pore_water_column(c)}_mM")]
    out_cols = [c for c in df.columns if c.startswith(tuple(MATRIC_COLS + CONDUCT_COLS))] + pw_cols
    df[out_cols].round(6).to_csv(DATA_DIR / "matric_conductivity_processed.csv")
    plot(df)
    plt.show()
