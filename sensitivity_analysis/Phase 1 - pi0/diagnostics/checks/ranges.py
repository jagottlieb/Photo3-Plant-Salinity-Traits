"""Analysis-window ranges of plotted variables in a batch (for choosing ylims)."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
ts = REPO / "sensitivity_analysis" / "Phase 1 - pi0" / "results" / (sys.argv[1] if len(sys.argv) > 1 else "20261002_1614") / "timeseries"
for f in sorted(ts.glob("*.csv")):
    df = pd.read_csv(f)
    an = df[~df["is_burn_in"]]
    qs = np.array([json.loads(c) for c in an["qs"]])
    ps = np.array([json.loads(c) for c in an["psi_s"]])
    print(f.name, f"qs per comp max {qs.max():.4f} min {qs.min():.4f}", f"psi_s min {ps.min():.4f}")
