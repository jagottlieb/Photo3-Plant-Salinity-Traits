#!/usr/bin/env python
"""
Run the targeted model runs defined in a sensitivity analysis Phase config.

Each Phase folder (e.g. ``sensitivity_analysis/Phase 1 - pi0/``) holds a
``config.yaml`` that lists baseline parameters, scenarios and named runs.
Every execution of this script creates a new batch folder:

    <phase>/results/<batch_id>/
        config.yaml        frozen copy of the config plus batch metadata
        manifest.csv       one row per run x scenario
        timeseries/<run>__<scenario>.csv   every model output variable

Usage (from the repo root):
    python sensitivity_analysis/sensitivity_analysis.py "sensitivity_analysis/Phase 1 - pi0/config.yaml"
"""

import argparse
import copy
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
# Model modules live in the repo root, one level above this script.
sys.path.insert(0, str(REPO_ROOT))

import hydraulics
import species_traits
from defs import Atmosphere, SimulationMultiComp, steps
from photosynthesis import C3_c_leaf_reduc
import soil
from soil import ConstantSoil, DrydownSoil, SaltySoilMultiple

# Scenario keys that are not model parameters (everything else must exist in baseline).
SCENARIO_ONLY_KEYS = {"label", "post_burn_cs", "post_burn_soil_dynamics"}
# Run keys that are metadata, not model parameters. ``group`` splits figures by run group.
RUN_META_KEYS = {"name", "group"}


# ---------------------------------------------------------------------------
# Config and batch bookkeeping
# ---------------------------------------------------------------------------

def load_config(path: Path) -> Dict[str, Any]:
    """Load a Phase config and check that scenario and run keys exist in baseline.

    Args:
        path: Path to the Phase ``config.yaml``.

    Returns:
        Parsed config dictionary.
    """
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    for section in ("simulation", "baseline", "scenarios", "active_scenarios", "runs"):
        if section not in cfg:
            raise KeyError(f"{path}: missing required section '{section}'")

    baseline_keys = set(cfg["baseline"])
    for name, scenario in cfg["scenarios"].items():
        unknown = set(scenario) - baseline_keys - SCENARIO_ONLY_KEYS
        if unknown:
            raise KeyError(f"Scenario {name}: keys not in baseline: {sorted(unknown)}")
    for name in cfg["active_scenarios"]:
        if name not in cfg["scenarios"]:
            raise KeyError(f"active_scenarios lists unknown scenario '{name}'")

    run_names = [run.get("name") for run in cfg["runs"]]
    if None in run_names or len(set(run_names)) != len(run_names):
        raise ValueError("Every run needs a unique 'name'")
    for run in cfg["runs"]:
        unknown = set(run) - RUN_META_KEYS - baseline_keys
        if unknown:
            raise KeyError(f"Run {run['name']}: keys not in baseline: {sorted(unknown)}")

    return cfg


def resolve_params(cfg: Dict[str, Any], run: Dict[str, Any], scenario_name: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Build the full parameter set for one run: baseline, then scenario overrides, then run overrides.

    Args:
        cfg: Parsed Phase config.
        run: One entry from ``cfg["runs"]``.
        scenario_name: Key into ``cfg["scenarios"]``.

    Returns:
        ``(params, scenario)``: model parameters, and the scenario dict
        (label and post-burn soil settings).
    """
    scenario = cfg["scenarios"][scenario_name]
    params = copy.deepcopy(cfg["baseline"])
    params.update({k: v for k, v in scenario.items() if k not in SCENARIO_ONLY_KEYS})
    params.update({k: v for k, v in run.items() if k not in RUN_META_KEYS})
    return params, scenario


def _git_info() -> Dict[str, Any]:
    """Return the current git commit hash and whether the working tree is dirty."""
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=REPO_ROOT, capture_output=True, text=True, check=True,
        ).stdout.strip())
        return {"git_commit": commit, "git_dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"git_commit": None, "git_dirty": None}


def create_batch(config_path: Path) -> Path:
    """Create ``results/<batch_id>/`` and write the frozen config for this batch.

    The frozen config is the original config text (comments kept) followed by a
    ``batch`` block with metadata only: batch id, timestamp and git state.

    Args:
        config_path: Path to the Phase ``config.yaml``.

    Returns:
        Path to the new batch folder.
    """
    created = datetime.now()
    batch_id = created.strftime("%Y%m%d_%H%M")
    batch_dir = config_path.parent / "results" / batch_id
    if batch_dir.exists():
        raise FileExistsError(f"Batch folder already exists: {batch_dir}")
    (batch_dir / "timeseries").mkdir(parents=True)

    batch_meta = {
        "batch": {
            "batch_id": batch_id,
            "created": created.isoformat(timespec="seconds"),
            **_git_info(),
        }
    }
    original_text = config_path.read_text(encoding="utf-8")
    with open(batch_dir / "config.yaml", "w", encoding="utf-8") as f:
        f.write(original_text.rstrip() + "\n\n")
        f.write("# --- Batch metadata (written by sensitivity_analysis.py) ---\n")
        yaml.safe_dump(batch_meta, f, sort_keys=False, default_flow_style=None)
    return batch_dir


# ---------------------------------------------------------------------------
# Model setup and execution
# ---------------------------------------------------------------------------

def load_weather(path: Path, timestepM: float, total_days: float) -> Dict[str, List[float]]:
    """Load weather forcing and linearly resample it to the model timestep.

    Args:
        path: Excel file with Temperature (C), Relative Humidity (%) and GHI (W/m2).
        timestepM: Model timestep (min).
        total_days: Burn-in plus analysis window (days); must fit in the file.

    Returns:
        Dict with ``phi`` (W/m2), ``ta`` (K) and ``qa`` (kg/kg) lists.
    """
    df = pd.read_excel(path, engine="openpyxl")
    required_cols = ["Temperature", "Relative Humidity", "GHI"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Weather file missing required columns: {missing_cols}")

    # Infer source timestep from a datetime column when available; otherwise assume 30 min.
    source_timestepM = 30.0
    for dt_col in ["DateTime", "Datetime", "datetime", "Timestamp", "timestamp", "Time"]:
        if dt_col in df.columns:
            t = pd.to_datetime(df[dt_col], errors="coerce").dropna()
            if len(t) > 1:
                delta_min = np.diff(t.values).astype("timedelta64[s]").astype(float) / 60.0
                delta_min = delta_min[np.isfinite(delta_min) & (delta_min > 0)]
                if len(delta_min) > 0:
                    source_timestepM = float(np.median(delta_min))
            break

    src_minutes = np.arange(len(df), dtype=float) * source_timestepM
    target_minutes = np.arange(int(np.floor(src_minutes[-1] / timestepM)) + 1, dtype=float) * timestepM

    needed = int(steps(total_days, int(timestepM)))
    if len(target_minutes) < needed:
        available_days = src_minutes[-1] / (24 * 60)
        raise ValueError(
            f"Weather covers {available_days:.1f} days but burn-in + analysis needs {total_days} days"
        )

    def resample(col: str) -> np.ndarray:
        return np.interp(target_minutes, src_minutes, df[col].to_numpy(dtype=float))

    temp_c = resample("Temperature")
    rh = resample("Relative Humidity")
    psat = 611.2 * np.exp((17.67 * temp_c) / (243.5 + temp_c))
    return {
        "phi": list(resample("GHI")),
        "ta": list(temp_c + 273.0),
        "qa": list(0.622 * (rh / 100.0) * psat / 101325),
    }


def resolve_zr_arr(sim_cfg: Dict[str, Any], params: Dict[str, Any], species_obj: Any) -> np.ndarray:
    """Return compartment depths from the config, or from the species rooting depth.

    ``simulation.zr_arr: species`` sets every compartment depth to ``species.ZR``,
    with one compartment per ``s_init`` entry. An explicit list is used as given.

    Args:
        sim_cfg: The config ``simulation`` section.
        params: Resolved model parameters (``s_init`` and ``root_frac`` set the compartment count).
        species_obj: Species instance from ``species_traits``.

    Returns:
        Compartment depths (m), one per compartment.
    """
    n_comp = len(params["s_init"])
    zr_cfg = sim_cfg["zr_arr"]
    if isinstance(zr_cfg, str):
        if zr_cfg.strip().lower() != "species":
            raise ValueError(f"zr_arr must be a list or 'species', got '{zr_cfg}'")
        zr_arr = np.full(n_comp, float(species_obj.ZR))
    else:
        zr_arr = np.array(zr_cfg, dtype=float)

    if len(zr_arr) != n_comp or len(params["root_frac"]) != n_comp:
        raise ValueError(
            f"zr_arr ({len(zr_arr)}), root_frac ({len(params['root_frac'])}) and "
            f"s_init ({n_comp}) must have the same number of compartments"
        )
    return zr_arr


def build_model(params: Dict[str, Any], sim_cfg: Dict[str, Any], weather: Dict[str, List[float]]) -> Dict[str, Any]:
    """Build a SimulationMultiComp model with stem and leaf storage from config values.

    Args:
        params: Resolved model parameters (baseline + scenario + run).
        sim_cfg: The config ``simulation`` section.
        weather: Output of :func:`load_weather` (first step initialises the atmosphere).

    Returns:
        Dict with ``plant``, ``hydro`` and ``species`` objects.
    """
    dt = float(sim_cfg["timestepM"]) * 60.0
    root_frac_arr = np.array(params["root_frac"], dtype=float)
    cs_init = np.array(params["cs_init"], dtype=float)

    species_obj = getattr(species_traits, sim_cfg["species"])()
    zr_arr = resolve_zr_arr(sim_cfg, params, species_obj)
    species_obj.GPMAX = species_obj.GPMAX * float(params.get("gpmax_scale", 1.0))
    species_obj.VWT = float(params["VWT"])
    species_obj.VWTLEAF = float(params["VWTLEAF"])
    species_obj.GWMAXLEAF = float(params["GWMAXLEAF"])

    gw_mode = str(sim_cfg.get("gwmax_mode", "equation")).strip().lower()
    if gw_mode == "manual":
        if params.get("GWMAX_manual") is None:
            raise ValueError("gwmax_mode is manual but baseline.GWMAX_manual is null")
        species_obj.GWMAX = float(params["GWMAX_manual"])
    elif gw_mode == "equation":
        species_obj.GWMAX = species_obj.CAP * species_obj.VWT * 10**6 / (0.63 * 4 * 60 * 60)
    else:
        raise ValueError(f"Unknown gwmax_mode: {gw_mode}")

    gcut = params.get("gcut", "species")
    gcut = float(species_obj.GCUT) if str(gcut).strip().lower() == "species" else float(gcut)

    atmosphere_obj = Atmosphere(weather["phi"][0], weather["ta"][0], weather["qa"][0])
    soil_obj = SaltySoilMultiple(
        stype=getattr(soil, sim_cfg.get("soil_texture", "Berger"))(),
        dynamics=DrydownSoil(),
        zr=zr_arr,
        s=np.array(params["s_init"], dtype=float),
        cs=cs_init.copy(),
    )
    photo_obj = C3_c_leaf_reduc(
        species_obj,
        atmosphere_obj,
        pi0_leaf=params["pi0_leaf"],
        eta_leaf=params["eta_leaf"],
    )
    hydro_obj = hydraulics.HalophyteStemLeafStorageMultiComp(
        species=species_obj,
        atm=atmosphere_obj,
        soil=soil_obj,
        photo=photo_obj,
        vwi_stem=params["vwi_stem"],
        c_stem=params["c_stem"],
        s_arr=np.array(params["s_init"], dtype=float),
        root_frac_arr=root_frac_arr,
        B=params["B_param"],
        cs_arr=cs_init.copy(),
        wr_stem=params["wr_stem"],
        wft_stem=params["wft_stem"],
        pi0_stem=params["pi0_stem"],
        eta_stem=params["eta_stem"],
        mcap_stem=params["mcap_stem"],
        dt=dt,
        salt_uptake=bool(params["salt_uptake"]),
        psi_wf_mode=sim_cfg.get("psi_wf_mode", "bartlett"),
        vwi_leaf=params["vwi_leaf"],
        c_leaf=params["c_leaf"],
        wr_leaf=params["wr_leaf"],
        wft_leaf=params["wft_leaf"],
        pi0_leaf=params["pi0_leaf"],
        eta_leaf=params["eta_leaf"],
        mcap_leaf=params["mcap_leaf"],
        leaf_uptake_frac=params["leaf_uptake_frac"],
        gcut=gcut,
        E=float(params["E"]),
        F_CAP=params["F_CAP"],
        dynamic_E=bool(params["dynamic_E"]),
        c_stem_max=params["c_stem_max"],
        kr=float(params.get("kr", 1e-8)),
    )
    photo_obj.hydro_ref = hydro_obj

    plant_obj = SimulationMultiComp(
        species_cls=species_obj,
        atm_cls=atmosphere_obj,
        soil_cls=soil_obj,
        photo_cls=photo_obj,
        hydro_cls=hydro_obj,
        zr_arr=zr_arr,
        root_frac_arr=root_frac_arr,
        B=params["B_param"],
        dt=dt,
    )
    return {"plant": plant_obj, "hydro": hydro_obj, "species": species_obj}


def run_simulation(
    model: Dict[str, Any],
    params: Dict[str, Any],
    scenario: Dict[str, Any],
    weather: Dict[str, List[float]],
    sim_cfg: Dict[str, Any],
    tag: str = "",
) -> Dict[str, Any]:
    """Step the model through burn-in and the analysis window.

    During burn-in soil moisture is held constant and salt uptake is off. At the
    end of burn-in the scenario's per-compartment soil dynamics and salt uptake
    setting take effect, and soil salinity is reset to ``post_burn_cs`` if enabled.

    Args:
        model: Output of :func:`build_model`.
        params: Resolved model parameters.
        scenario: Scenario dict (``post_burn_cs``, ``post_burn_soil_dynamics``).
        weather: Output of :func:`load_weather`.
        sim_cfg: The config ``simulation`` section.
        tag: Label used in progress messages.

    Returns:
        Dict with ``output`` (``plant.output()``), ``n_steps`` and ``burn_in_steps``.
    """
    plant = model["plant"]
    hydro = model["hydro"]
    timestepM = int(sim_cfg["timestepM"])
    dt = timestepM * 60.0
    burn_in_steps = int(steps(float(sim_cfg["burn_in_days"]), timestepM))
    n_steps = int(steps(float(sim_cfg["burn_in_days"]) + float(sim_cfg["analysis_days"]), timestepM))
    n_comp = len(plant.soil.s)

    post_dynamics = [
        DrydownSoil() if str(mode).strip().lower() == "drydown" else ConstantSoil()
        for mode in scenario.get("post_burn_soil_dynamics", ["constant"] * n_comp)
    ]
    post_cs = np.array(scenario.get("post_burn_cs", params["cs_init"]), dtype=float)
    if len(post_dynamics) != n_comp or len(post_cs) != n_comp:
        raise ValueError(f"{tag}: post_burn settings must have {n_comp} entries")

    plant.soil.dynamics = [ConstantSoil() for _ in range(n_comp)]
    hydro.Salt_Uptake = False
    progress_stride = max(1, n_steps // 4)

    for i in range(n_steps):
        if i == burn_in_steps:
            plant.soil.dynamics = post_dynamics
            hydro.Salt_Uptake = bool(params["salt_uptake"])
            if sim_cfg.get("enable_post_burn_cs_reset", True):
                plant.soil.cs = post_cs.copy()
                if hasattr(plant.soil, "MS"):
                    plant.soil.MS = plant.soil.cs * plant.soil.ZR * plant.soil.N * plant.soil.s
        plant.update(dt, weather["phi"][i], weather["ta"][i], weather["qa"][i])
        if (i + 1) % progress_stride == 0 or i == n_steps - 1:
            print(f"  {tag}: {i + 1}/{n_steps} steps ({(i + 1) / n_steps * 100:.0f}%)")

    return {"output": plant.output(), "n_steps": n_steps, "burn_in_steps": burn_in_steps}


def outputs_to_dataframe(result: Dict[str, Any], weather: Dict[str, List[float]], sim_cfg: Dict[str, Any]) -> pd.DataFrame:
    """Put every model output variable into one row per timestep.

    Per-compartment variables (e.g. ``qs``, ``s``) keep the model's pair form:
    each cell holds a JSON list with one value per compartment. Scalar outputs
    become constant columns. Compatibility aliases that point to the same list
    as an earlier key (e.g. ``qw`` for ``qw_stem``) are skipped. Adds the
    weather forcing and time columns.

    Args:
        result: Output of :func:`run_simulation`.
        weather: Output of :func:`load_weather`.
        sim_cfg: The config ``simulation`` section.

    Returns:
        DataFrame with ``n_steps`` rows.
    """
    n = result["n_steps"]
    burn_in_steps = result["burn_in_steps"]
    dt = float(sim_cfg["timestepM"]) * 60.0
    step = np.arange(n)
    time_days = (step + 1) * dt / 86400.0

    columns: Dict[str, Any] = {
        "step": step,
        "time_days": time_days,
        "days_since_burn_in": time_days - burn_in_steps * dt / 86400.0,
        "is_burn_in": step < burn_in_steps,
        "forcing_phi": weather["phi"][:n],
        "forcing_ta": weather["ta"][:n],
        "forcing_qa": weather["qa"][:n],
    }

    seen_ids = set()
    for key, value in result["output"].items():
        if id(value) in seen_ids:
            continue
        seen_ids.add(id(value))
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0:
            columns[key] = np.full(n, float(arr))
        elif arr.ndim == 1 and arr.shape[0] == n:
            columns[key] = arr
        elif arr.ndim == 2 and arr.shape[1] == n:
            columns[key] = [json.dumps(row) for row in arr.T.tolist()]
        elif arr.ndim == 2 and arr.shape[0] == n:
            columns[key] = [json.dumps(row) for row in arr.tolist()]
        else:
            print(f"  Skipping output '{key}' with unexpected shape {arr.shape}")

    return pd.DataFrame(columns)


# ---------------------------------------------------------------------------
# Command line entry point
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> Path:
    """Run every run x active scenario in a Phase config and save the results.

    Returns:
        Path to the batch folder.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="Path to the Phase config.yaml")
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    cfg = load_config(config_path)

    sim_cfg = cfg["simulation"]
    total_days = float(sim_cfg["burn_in_days"]) + float(sim_cfg["analysis_days"])
    weather = load_weather(REPO_ROOT / sim_cfg["weather_file"], float(sim_cfg["timestepM"]), total_days)

    batch_dir = create_batch(config_path)
    print(f"Batch folder: {batch_dir}")

    manifest_rows = []
    n_total = len(cfg["active_scenarios"]) * len(cfg["runs"])
    k = 0
    for scenario_name in cfg["active_scenarios"]:
        for run in cfg["runs"]:
            k += 1
            tag = f"{run['name']}__{scenario_name}"
            print(f"\n[{k}/{n_total}] {tag}")
            params, scenario = resolve_params(cfg, run, scenario_name)
            start = time.perf_counter()
            model = build_model(params, sim_cfg, weather)
            result = run_simulation(model, params, scenario, weather, sim_cfg, tag=tag)
            df = outputs_to_dataframe(result, weather, sim_cfg)
            csv_name = f"timeseries/{tag}.csv"
            df.to_csv(batch_dir / csv_name, index=False)
            runtime_s = time.perf_counter() - start

            manifest_rows.append({
                "run": run["name"],
                "group": run.get("group"),
                "overrides": json.dumps({k_: v for k_, v in run.items() if k_ not in RUN_META_KEYS}),
                "scenario": scenario_name,
                "scenario_label": scenario.get("label", scenario_name),
                "n_steps": result["n_steps"],
                "burn_in_steps": result["burn_in_steps"],
                "runtime_s": round(runtime_s, 1),
                "timeseries_file": csv_name,
            })
            pd.DataFrame(manifest_rows).to_csv(batch_dir / "manifest.csv", index=False)
            print(f"  Saved {csv_name} ({runtime_s:.0f} s)")

    print(f"\nDone: {n_total} runs in {batch_dir}")
    return batch_dir


if __name__ == "__main__":
    main()
