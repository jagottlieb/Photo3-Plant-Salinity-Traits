"""s_init sweep, round 2: run Phase 1 S2 (drydown) for every pi0 run at several s_init values.

Ran with 0.2 m compartments (zr_arr: species), before the gp lag fix in hydraulics.py,
so E is somewhat inflated (most for the high run). Values used: 0.34, 0.36, 0.38, then
0.40, 0.42, 0.45 (pass them as arguments). Writes one timeseries CSV per (s_init, run)
next to this script without creating a batch folder. Config is read but never modified.
"""

import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents if (p / "hydraulics.py").exists())
sys.path.insert(0, str(REPO / "SensitivityAnalysis"))

from sensitivity_analysis import (  # noqa: E402
    build_model, load_config, load_weather, outputs_to_dataframe, resolve_params, run_simulation,
)

CONFIG = REPO / "SensitivityAnalysis" / "Phase 1 - pi0" / "config.yaml"
OUT = Path(__file__).resolve().parent
S_INITS = [0.34, 0.36, 0.38]


def one(args):
    s_init, run_name = args
    cfg = load_config(CONFIG)
    sim = cfg["simulation"]
    weather = load_weather(REPO / sim["weather_file"], float(sim["timestepM"]),
                           float(sim["burn_in_days"]) + float(sim["analysis_days"]))
    run = next(r for r in cfg["runs"] if r["name"] == run_name)
    params, scenario = resolve_params(cfg, run, "S2")
    params["s_init"] = [s_init, s_init]
    t0 = time.perf_counter()
    model = build_model(params, sim, weather)
    result = run_simulation(model, params, scenario, weather, sim, tag=f"s{s_init}_{run_name}")
    df = outputs_to_dataframe(result, weather, sim)
    df.to_csv(OUT / f"s{int(round(s_init * 100))}_{run_name}.csv", index=False)
    return s_init, run_name, time.perf_counter() - t0


if __name__ == "__main__":
    s_inits = [float(a) for a in sys.argv[1:]] or S_INITS
    jobs = [(s, r) for s in s_inits for r in ("low", "med", "high")]
    with ProcessPoolExecutor(max_workers=5) as ex:
        for s, r, dt in ex.map(one, jobs):
            print(f"done s_init={s} {r} ({dt:.0f} s)")
