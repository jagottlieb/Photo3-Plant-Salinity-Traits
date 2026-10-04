# Systematic Sensitivity Analysis Plan

## Overview
The initial sensitivity analysis for the Photo3-Plant-Salinity-Traits model is organised as a series of **phases**. Each phase is a targeted set of model runs (one or more parameters varied, everything else held at baseline), defined by its own config file. Every execution of a phase is saved as a batch with a frozen copy of the config, so there is a record of exactly what was tested.

---

## 1. Parameter Categories

### Stem Storage Parameters (9)
- `pi0_stem`: Stem osmotic potential at full turgor (MPa)
- `wft_stem`: Relative water content at full turgor
- `wr_stem`: Residual water fraction (relative to VWT)
- `eta_stem`: Bulk elastic modulus (MPa)
- `mcap_stem`: Slope of psi_w above full turgor (MPa)
- `vwi_stem`: Initial water storage fraction
- `c_stem`: Initial stem storage salt concentration (mol/m3)
- `VWT`: Stem storage volume (m³/m² leaf)
- `GWMAX`: Stem storage conductance (um/MPa/s)

### Leaf Storage Parameters (9)
- `pi0_leaf`: Leaf osmotic potential at full turgor (MPa)
- `wft_leaf`: Relative water content at full turgor
- `wr_leaf`: Residual water fraction (relative to VWTLEAF)
- `eta_leaf`: Bulk elastic modulus (MPa)
- `mcap_leaf`: Slope of psi_w above full turgor (MPa)
- `vwi_leaf`: Initial water storage fraction
- `c_leaf`: Initial leaf storage salt concentration (mol/m3)
- `VWTLEAF`: Leaf storage volume (m³/m² leaf)
- `GWMAXLEAF`: Leaf storage conductance (um/MPa/s)

### Soil Parameters (5)
- `s_init`: Initial relative soil moisture per compartment
- `cs_init`: Initial soil salinity per compartment (mM)
- `root_frac`: Root fraction per compartment
- `B_param`: Root length density parameter
- `F_CAP`: Relative position of the stem storage node along the hydraulic path

### Salt and Conductance Parameters (4)
- `E`: Salt filtration efficiency
- `dynamic_E`: Dynamic filtration efficiency
- `salt_uptake`: Salt uptake enabled
- `leaf_uptake_frac`: Leaf share of salt uptake

---

## 2. Target Outputs

Every batch saves **all** model output variables per timestep (see `timeseries/*.csv`). The ones of main interest are:

### Water Status
- `psi_l`: Leaf water potential (MPa)
- `psi_w_stem`, `psi_w_leaf`: Storage water potentials (MPa)
- `psi_w_osm_leaf`, `psi_w_turgor_leaf`: Leaf osmotic and turgor components (MPa)

### Water Fluxes (um/s, ground area basis)
- `ev`: Transpiration
- `qs`: Root water uptake, saved as a per-compartment pair (as are `s`, `psi_s`, `cs`, `gsr`)
- `qw_stem`, `qw_leaf`: Storage-to-xylem fluxes (positive = release to xylem)

### Photosynthesis
- `a`: Photosynthesis rate (umol/m²/s)
- `gsw`: Stomatal conductance

---

## 3. Phases

| Phase | Folder | Varied | Runs |
|-------|--------|--------|------|
| 1 | `SensitivityAnalysis/Phase 1 - pi0/` | `pi0_leaf` and `pi0_stem`, covaried low/med/high | 3 per scenario |
| 2+ | to be defined | | |

### Phase 1 - pi0
- Low: `pi0_leaf` = -1.1, `pi0_stem` = -0.9
- Med (baseline): `pi0_leaf` = -1.4, `pi0_stem` = -1.1
- High: `pi0_leaf` = -1.7, `pi0_stem` = -1.3
- Scenario: Drydown (S2) by default; set `active_scenarios` in the config. S2 starts both compartments at `s_init` = 0.42 and dries both after burn-in. Compartment depths are the species rooting depth (`zr_arr: species`, 0.2 m for Pecan), so soil-root conductance and the soil water balance use the same depth. With 0.2 m compartments s falls from 0.42 to about 0.255 over 15 days; transpiration declines from about day 3 and reaches its cuticular floor on days 11-13, so the pi0 runs stay separate for most of the window (starting from 0.32 they merged by day 8; s_init sweep 0.32-0.45). (Earlier 1 m runs: starting from 0.45 barely dried the soil, and starting below about 0.29 put transpiration at its floor from day 0; starting below 0.26 made the leaf water potential solver fail)
- Simulation: 7-day burn-in plus 15-day analysis window (weather is `Lab_Weather_Summer_Average_Day_30min.xlsx`: the time-of-day mean of the lab record, repeated for 45 days)
- Figures: 15-day time series of transpiration, root uptake, soil water potential, storage fluxes, psi_l and flux partition (share of E from roots, stem and leaf storage)
- Final batch: `20261002_1614` (after the gp lag fix in `hydraulics.py`)

#### Phase 1 diagnostics
Scripts used to choose the S2 drydown and explain the model behaviour. They are stowed for the record, not part of the standard workflow. Each script's docstring states the conditions it ran under; outputs are written next to the script. Both s_init sweeps and the diagnostics below ran before the gp lag fix, so their E values are somewhat inflated (most for the high run).
- `diagnostics/s_init_sweep/round1_zr1m/`: 1 m compartments. `drydown_sweep.py` (med run, s_init 0.30-0.23) and `drydown_pick.py` (all runs, 0.32/0.30/0.29), which led to s_init 0.32.
- `diagnostics/s_init_sweep/round2_zr02m/`: 0.2 m compartments. `run_sweep.py` (all runs, s_init 0.34-0.45; CSVs gitignored) and `summarize_sweep.py` (table and comparison figure), which led to s_init 0.42.
- `diagnostics/root_conductance/`: soil-root conductance vs soil water, the psi_l collapse at s_init 0.23, and the s = 0.45 vs 0.32 day-1 comparison.
- `diagnostics/night_psi_l/`: why night psi_l stays far below psi_s once the soil dries (stem storage sets it).
- `diagnostics/checks/`: quick checks (compartment depths, batch summary, plot ranges).
- `extra_plots/flux_partition/`: daily supply shares and prototype flux-partition figures (the plot itself is now `flux_partition` in `post_processing.py`).

Candidate parameters for later phases: storage volumes (`VWT`, `VWTLEAF`), storage conductances (`GWMAX`, `GWMAXLEAF`), elastic moduli (`eta_stem`, `eta_leaf`), and salinity settings (`cs_init`, `E`, `dynamic_E`).

---

## 4. Workflow

```
python SensitivityAnalysis/sensitivity_analysis.py "SensitivityAnalysis/Phase 1 - pi0/config.yaml"
python SensitivityAnalysis/post_processing.py "SensitivityAnalysis/Phase 1 - pi0/config.yaml" --batch latest
```

- `sensitivity_analysis.py` runs every run x active scenario and saves all outputs. Everything it does is set in the config.
- `post_processing.py` has one function per figure type. It reads the `plots` section of the **current** config, so axis limits and labels can be changed and figures re-made without re-running the model.

### Phase config sections
- `simulation`: species, weather file, timestep, burn-in and analysis days, compartment depths (`zr_arr`: a list in m, or `species` for the species rooting depth in every compartment), `gwmax_mode`
- `baseline`: every model parameter, with source comments
- `scenarios` and `active_scenarios`
- `runs`: named overrides of baseline values (add rows for more combinations)
- `plots`: shared settings plus one block per figure

### Folder layout
```
SensitivityAnalysis/
  Phase 1 - pi0/
    config.yaml
    results/
      <batch_id>/
        config.yaml         frozen config + batch metadata (id, time, git commit) (tracked in git)
        manifest.csv        one row per run x scenario, with each run's overrides
        timeseries/         <run>__<scenario>.csv, all output variables
        figures/<scenario>/ transpiration, root_uptake, storage_fluxes, psi_l
    diagnostics/            optional: scripts used to tune or explain the phase (tracked)
    extra_plots/            optional: scripts for figures beyond post_processing.py (tracked)
```
Only the frozen `config.yaml` in each batch is tracked; results and figures are gitignored.

To start a new phase, copy a Phase folder's `config.yaml` into a new `Phase N - <name>/` folder and edit the `runs` list.

---

## 5. Notes

- Each run takes about 8 s (22 days at 30 min steps), so a phase finishes in under a minute per scenario.
- Storage volumes and conductances match the `Pecan` class in `species_traits.py` (`VWTLEAF` = 5e-5 and `GWMAXLEAF` = 0.0005, both from American beech). Other baseline values match `analysis_notebook_sensitivity_v3.ipynb`, including `eta_leaf` = 17.5 and `E` = 0.95 (the `hydraulics.py` defaults of 5 and 0.99 are not used).

---

*Created: 2026-09-29. Restructured into phases: 2026-10-01*
