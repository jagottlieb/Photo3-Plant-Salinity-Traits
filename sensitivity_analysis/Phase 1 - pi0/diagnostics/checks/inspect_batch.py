"""Print S2 summary per run for a Phase 1 batch: psi_l at the -10 bound, soil drydown, daily peak E and min psi_l.

Usage: python inspect_batch.py <batch_id>
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PHASE = Path(__file__).resolve().parents[2]
batch = PHASE / "results" / sys.argv[1]
for run in ["low", "med", "high"]:
    df = pd.read_csv(batch / "timeseries" / f"{run}__S2.csv")
    a = df[~df.is_burn_in].copy()
    s = np.array([json.loads(x) for x in a["s"]])
    ps = np.array([json.loads(x) for x in a["psi_s"]])
    qs = np.array([json.loads(x) for x in a["qs"]])
    d = a["days_since_burn_in"].to_numpy()
    day = np.floor(d - 1e-9).astype(int)
    pinned = a["psi_l"] <= -9.99
    print(f"== {run}")
    print(f"  psi_l min {a.psi_l.min():.3f}; steps at -10: {pinned.sum()}"
          + (f", first at day {d[pinned.to_numpy()][0]:.2f}, s={s[pinned.to_numpy()][0]}" if pinned.any() else ""))
    print(f"  s start {s[0]}, end {s[-1]}; psi_s end {ps[-1]}")
    print(f"  qs max {qs.max(axis=0)}; psi_s min {ps.min(axis=0)}")
    for k in [0, 1, 2, 3, 5, 7, 10, 14]:
        m = day == k
        if m.any():
            print(f"  day {k+1:2d}: peak E {a.ev[m].max():.4f}, min psi_l {a.psi_l[m].min():.3f}, s_end {s[m][-1].round(4)}, psi_s {ps[m][-1].round(4)}")
