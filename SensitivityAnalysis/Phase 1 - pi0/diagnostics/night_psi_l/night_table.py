"""Pre-dawn and daytime path diagnostics for one S2 timeseries (night psi_l question).

Defaults to batch 20261002_1550 (s_init 0.32, 0.2 m compartments, before the gp lag fix).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
BATCH = REPO / "SensitivityAnalysis" / "Phase 1 - pi0" / "results" / "20261002_1550"
LAI, F_CAP = 17.0, 0.5

path = Path(sys.argv[1]) if len(sys.argv) > 1 else BATCH / "timeseries" / "med__S2.csv"
df = pd.read_csv(path)
pair = lambda col: np.array([json.loads(c) for c in df[col]])
qs, gsr, psi_s, s = pair("qs"), pair("gsr"), pair("psi_s"), pair("s")
df["qs_sum"], df["gsr_tot"], df["psi_s1"], df["s1"] = qs.sum(1), gsr.sum(1), psi_s[:, 0], s[:, 0]
df["gpL"] = df["gp"] * LAI
df["hour"] = np.round((df["time_days"] % 1) * 24, 2)
df["day"] = np.floor(df["days_since_burn_in"] + 1e-9).astype(int)

# gradients along the path
df["d_soil"] = df["psi_s1"] - df["psi_b"]                     # soil -> root base, across gsr
df["d_bx"] = df["psi_b"] - df["psi_x"]                        # root base -> stem node, across 2 gp LAI
df["d_xl"] = df["psi_x"] - df["psi_l"]                        # stem node -> leaf, across 2 gp LAI
df["d_wx"] = df["psi_w_stem"] - df["psi_x"]                   # stem storage -> stem node, across gw LAI
df["qs_over_gsr"] = df["qs_sum"] / df["gsr_tot"]
df["qbx_minus_qs"] = df["qbx"] - df["qs_sum"]
df["gp_lag_check"] = df["qbx"] + df["qw_stem"] + df["qw_leaf"] - df["ev"]
cols = ["day", "hour", "s1", "psi_s1", "psi_b", "psi_x", "psi_l", "psi_w_stem", "psi_w_leaf", "w_stem",
        "ev", "qs_sum", "qbx", "qw_stem", "qw_leaf", "gsr_tot", "gpL", "gw_stem",
        "d_soil", "d_bx", "d_xl", "d_wx", "qs_over_gsr", "flux_balance"]
pd.set_option("display.width", 400, "display.max_columns", 50, "display.precision", 4)
an = df[df["days_since_burn_in"] > -1]
print("PRE-DAWN (05:00)")
print(an[an["hour"] == 5.0][cols].to_string(index=False))
print("\nMIDDAY (13:00)")
print(an[an["hour"] == 13.0][cols].to_string(index=False))
print("\nflux_balance abs max/median:", df["flux_balance"].abs().max(), df["flux_balance"].abs().median())
print("qbx - sum qs abs max:", df["qbx_minus_qs"].abs().max())
print("(qbx+qw-ev) abs max:", df["gp_lag_check"].abs().max())
print("gw_stem*LAI range:", (df["gw_stem"] * LAI).min(), (df["gw_stem"] * LAI).max())
