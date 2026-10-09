#!/usr/bin/env python
"""
Run the model against a lab experiment window defined in a config.

Initial soil moisture and salinity come from the start-date probe data, and
irrigation and salt additions are applied as pulses (see ``lab_inputs.py``).
Model building, weather loading and batch bookkeeping are shared with
``sensitivity_analysis/sensitivity_analysis.py``. Every execution creates:

    results/<batch_id>/
        config.yaml               frozen copy of the config plus batch metadata
        timeseries/<name>.csv     every model output variable, with a datetime column

Usage (from the repo root):
    python lab_experiment_simulation/run_simulation.py lab_experiment_simulation/config.yaml
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "sensitivity_analysis"))

import sensitivity_analysis as sa
import soil
from defs import steps
from soil import ConstantSoil, DrydownSoil

import lab_inputs


def model_times(sim_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Start/end dates, step counts and step start times for the burn-in plus comparison window."""
    timestepM = int(sim_cfg["timestepM"])
    start = pd.Timestamp(sim_cfg["start_date"])
    end = pd.Timestamp(sim_cfg["end_date"]) + pd.Timedelta(days=1)
    burn_in_days = float(sim_cfg["burn_in_days"])
    total_days = burn_in_days + (end - start) / pd.Timedelta(days=1)
    n_steps = int(steps(total_days, timestepM))
    step_starts = start - pd.Timedelta(days=burn_in_days) + pd.to_timedelta(range(n_steps), unit="min") * timestepM
    return {
        "start": start,
        "end": end,
        "total_days": total_days,
        "n_steps": n_steps,
        "burn_in_steps": int(steps(burn_in_days, timestepM)),
        "step_starts": step_starts,
    }


def run_lab_simulation(model: Dict[str, Any], params: Dict[str, Any], weather: Dict[str, List[float]],
                       sim_cfg: Dict[str, Any], times: Dict[str, Any],
                       pulses: Dict[int, Dict[str, Any]]) -> Dict[str, Any]:
    """Step the model through burn-in and the comparison window, applying lab pulses.

    During burn-in soil moisture is held constant and salt uptake is off. From
    ``start_date`` the soil dries down, salt uptake follows ``params``, and each
    pulse is added at the start of its step.

    Returns:
        Dict with ``output`` (``plant.output()``), ``n_steps`` and ``burn_in_steps``.
    """
    plant, hydro = model["plant"], model["hydro"]
    dt = float(sim_cfg["timestepM"]) * 60.0
    n_steps, burn_in_steps = times["n_steps"], times["burn_in_steps"]
    n_comp = len(plant.soil.s)

    plant.soil.dynamics = [ConstantSoil() for _ in range(n_comp)]
    hydro.Salt_Uptake = False
    progress_stride = max(1, n_steps // 4)

    for i in range(n_steps):
        if i == burn_in_steps:
            plant.soil.dynamics = [DrydownSoil() for _ in range(n_comp)]
            hydro.Salt_Uptake = bool(params["salt_uptake"])
        if i in pulses:
            lab_inputs.apply_pulse(plant.soil, pulses[i]["ds"], pulses[i]["dms"])
        plant.update(dt, weather["phi"][i], weather["ta"][i], weather["qa"][i])
        if (i + 1) % progress_stride == 0 or i == n_steps - 1:
            print(f"  {i + 1}/{n_steps} steps ({(i + 1) / n_steps * 100:.0f}%)")

    return {"output": plant.output(), "n_steps": n_steps, "burn_in_steps": burn_in_steps}


def simulate(cfg: Dict[str, Any]) -> pd.DataFrame:
    """Set up lab inputs, run the model for a loaded config, and return the time series.

    Args:
        cfg: Parsed lab simulation config (``simulation``, ``lab`` and ``params`` sections).

    Returns:
        Output of ``sa.outputs_to_dataframe`` with a leading ``datetime`` column (step end times).
    """
    sim_cfg, lab, params = dict(cfg["simulation"]), cfg["lab"], dict(cfg["params"])

    times = model_times(sim_cfg)
    porosity = getattr(soil, sim_cfg["soil_texture"])().N
    zr = lab_inputs.compartment_depth(lab)
    sim_cfg["zr_arr"] = [zr] * len(lab["sides"])
    params["s_init"] = lab_inputs.initial_s(lab, times["start"], porosity)
    params["cs_init"] = lab_inputs.initial_cs(lab, times["start"])
    pulses = lab_inputs.pulse_schedule(lab, times["start"], times["end"], porosity, times["step_starts"])

    print(f"Compartment depth: {zr:.4f} m")
    print(f"Initial s:  {[round(v, 3) for v in params['s_init']]}")
    print(f"Initial cs: {[round(v, 3) for v in params['cs_init']]} mM")
    for step, p in sorted(pulses.items()):
        print(f"Pulse {p['time']} (step {step}): ds = {p['ds'].round(4).tolist()}, "
              f"salt = {p['dms'].round(4).tolist()} mol/m2")

    weather = sa.load_weather(REPO_ROOT / sim_cfg["weather_file"], float(sim_cfg["timestepM"]), times["total_days"])
    model = sa.build_model(params, sim_cfg, weather)
    result = run_lab_simulation(model, params, weather, sim_cfg, times, pulses)
    df = sa.outputs_to_dataframe(result, weather, sim_cfg)
    df.insert(0, "datetime", times["step_starts"] + pd.Timedelta(minutes=int(sim_cfg["timestepM"])))
    return df


def main(argv: Optional[List[str]] = None) -> Path:
    """Run the lab simulation in a config and save the results.

    Returns:
        Path to the batch folder.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="Path to the lab simulation config.yaml")
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    sim_cfg = cfg["simulation"]

    batch_dir = sa.create_batch(config_path)
    print(f"Batch folder: {batch_dir}")
    start = time.perf_counter()
    df = simulate(cfg)
    csv_path = batch_dir / "timeseries" / f"{sim_cfg['name']}.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved {csv_path.relative_to(batch_dir)} ({time.perf_counter() - start:.0f} s)")
    return batch_dir


if __name__ == "__main__":
    main()
