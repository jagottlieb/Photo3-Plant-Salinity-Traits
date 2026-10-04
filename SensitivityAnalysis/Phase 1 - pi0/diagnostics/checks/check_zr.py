"""Check that `zr_arr: species` gives 0.2 m (Pecan.ZR) in soil, plant and hydraulics, and that explicit lists still work."""

import sys
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
sys.path.insert(0, str(REPO / "SensitivityAnalysis"))
import sensitivity_analysis as sa

cfg = sa.load_config(REPO / "SensitivityAnalysis/Phase 1 - pi0/config.yaml")
sim = cfg["simulation"]
w = sa.load_weather(REPO / sim["weather_file"], float(sim["timestepM"]), 22)
p, sc = sa.resolve_params(cfg, cfg["runs"][1], "S2")
m = sa.build_model(p, sim, w)
print("soil.ZR", m["plant"].soil.ZR, "plant.zr_arr", m["plant"].zr_arr,
      "hydro.zr", m["hydro"].zr, "species.ZR", m["species"].ZR)
print("explicit", sa.resolve_zr_arr(dict(sim, zr_arr=[1.0, 1.0]), p, m["species"]))
