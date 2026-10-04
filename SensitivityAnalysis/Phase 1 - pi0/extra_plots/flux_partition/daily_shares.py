"""Print daily supply shares (% of daily E) for batch 20261002_1550 and the s_init=0.42 sweep runs.

Both sources predate the gp lag fix, so root + stem + leaf shares sum to less than 100%.
The figures here were made with the prototype of plot_flux_partition, now in post_processing.py.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

PHASE = Path(__file__).resolve().parents[2]
SRC = {"b1550": PHASE / "results" / "20261002_1550" / "timeseries" / "{r}__S2.csv",
       "s042": PHASE / "diagnostics" / "s_init_sweep" / "round2_zr02m" / "s42_{r}.csv"}
for tag, pat in SRC.items():
    for r in ("low", "med", "high"):
        df = pd.read_csv(str(pat).format(r=r))
        df = df[~df["is_burn_in"]]
        day = np.floor(df["days_since_burn_in"] - 1e-9).astype(int) + 1
        qs = df["qs"].map(lambda c: sum(json.loads(c)))
        g = pd.DataFrame({"day": day, "E": df["ev"], "root": qs, "stem": df["qw_stem"], "leaf": df["qw_leaf"]}).groupby("day").sum()
        sh = 100 * g[["root", "stem", "leaf"]].div(g["E"], axis=0)
        sh["sum"] = sh.sum(axis=1)
        print(tag, r, " | ".join(f"d{d}: root {sh.loc[d,'root']:.0f} stem {sh.loc[d,'stem']:.0f} sum {sh.loc[d,'sum']:.0f}" for d in (1, 5, 8, 10, 15)))
