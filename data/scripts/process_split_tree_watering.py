"""Convert the pecan Watering Chart into water and NaCl inputs for the split-root trees.

Reads the 'Watering Chart' sheet of Pecan_Data_Second_Round_May-July.xlsx and writes a
daily timeseries (one row per sheet date, stamped at WATERING_HOUR) with, for each split
tree and side, the irrigation volume (L) and moles of NaCl added to that bucket.

The 'Salt Conc Splits' value (dS/m) on a row is applied to every split bucket watered on
that row; rows without a value are treated as fresh water. Rows that look inconsistent
(salt listed but the note says fresh, or a note with no recorded split volumes) are
printed for manual review.
"""

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "tree_and_watering_data"
INPUT_FILE = DATA_DIR / "Pecan_Data_Second_Round_May-July.xlsx"
SHEET = "Watering Chart"
OUTPUT_FILE = DATA_DIR / "split_tree_water_salt_inputs.csv"

N_SPLIT_TREES = 4
SIDES = ("left", "right")
WATERING_HOUR = 14  # all watering events assumed to occur at 14:00
MOLES_DECIMALS = 4

# Electrical conductivity to salt concentration: 1 dS/m = 640 mg/L
# https://prometheusprotocols.net/experimental-design-and-analysis/experimental-treatments/salinity/
EC_TO_MG_PER_L = 640.0
# Molar mass of NaCl (g/mol), used to convert mg/L to mMol (mmol/L)
NACL_MOLAR_MASS = 58.44
# Volumes in the sheet are mL; 1 mL = 1e-3 L
ML_TO_L = 1e-3
# 1 mmol = 1e-3 mol
MMOL_TO_MOL = 1e-3


def ec_to_mmol(ec_ds_m: pd.Series) -> pd.Series:
    """Convert electrical conductivity (dS/m) to NaCl concentration (mMol = mmol/L)."""
    return ec_ds_m * EC_TO_MG_PER_L / NACL_MOLAR_MASS


def nacl_moles(vol_l: pd.Series, conc_mmol: pd.Series) -> pd.Series:
    """Moles of NaCl in vol_l liters of solution at conc_mmol (mmol/L)."""
    return conc_mmol * vol_l * MMOL_TO_MOL


def load_watering_chart(path: Path = INPUT_FILE) -> pd.DataFrame:
    """Load the Watering Chart with split-tree buckets named tree{i}_{side}."""
    raw = pd.read_excel(path, sheet_name=SHEET, header=None)
    header, sides = raw.iloc[0], raw.iloc[1]
    body = raw.iloc[2:].reset_index(drop=True)

    df = pd.DataFrame({
        "date": pd.to_datetime(body[0]),
        "notes": body[1],
        "salt_ec_splits": pd.to_numeric(body[2], errors="coerce"),
    })
    split_cols = [c for c in raw.columns if str(header[c]).startswith("Split Tree")]
    for c in split_cols:
        tree = int(str(header[c]).split()[-1])
        side = str(sides[c]).strip().lower()
        df[f"tree{tree}_{side}"] = pd.to_numeric(body[c], errors="coerce")
    return df


def flag_suspect_rows(df: pd.DataFrame, buckets: list[str]) -> None:
    """Print rows whose notes conflict with the recorded split-tree volumes or salt."""
    has_notes = df["notes"].notna()
    no_split_vol = df[buckets].isna().all(axis=1)
    fresh_note = df["notes"].astype(str).str.contains("fresh", case=False)
    has_salt = df["salt_ec_splits"].fillna(0) > 0

    checks = {
        "salt EC listed but note mentions fresh": has_salt & fresh_note,
        "note present but no split volumes recorded": has_notes & no_split_vol,
    }
    for label, mask in checks.items():
        if mask.any():
            print(f"Review ({label}):")
            print(df.loc[mask, ["date", "notes", "salt_ec_splits"]].to_string(index=False))


def build_inputs(df: pd.DataFrame) -> pd.DataFrame:
    """Return the 16-column timeseries of bucket volume (L) and NaCl (mol)."""
    buckets = [f"tree{t}_{s}" for t in range(1, N_SPLIT_TREES + 1) for s in SIDES]
    flag_suspect_rows(df, buckets)

    conc_mmol = ec_to_mmol(df["salt_ec_splits"].fillna(0.0))
    out = pd.DataFrame(index=df["date"] + pd.Timedelta(hours=WATERING_HOUR))
    out.index.name = "timestamp"
    for b in buckets:
        vol_l = df[b].fillna(0.0) * ML_TO_L
        out[f"{b}_vol_liters"] = vol_l.to_numpy()
        out[f"{b}_moles"] = nacl_moles(vol_l, conc_mmol).round(MOLES_DECIMALS).to_numpy()
    return out


def main() -> None:
    """Process the Watering Chart and write the split-tree input timeseries to CSV."""
    inputs = build_inputs(load_watering_chart())
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    inputs.to_csv(OUTPUT_FILE)
    print(f"Wrote {inputs.shape[0]} rows x {inputs.shape[1]} columns to {OUTPUT_FILE}")
    print(inputs[(inputs != 0).any(axis=1)].to_string(float_format="{:.4g}".format))


if __name__ == "__main__":
    main()
