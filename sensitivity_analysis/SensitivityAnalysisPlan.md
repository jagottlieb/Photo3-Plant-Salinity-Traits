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
- `leaf_uptake_frac`: Partition coefficient C_l, the leaf-to-stem ratio of salt concentration increases (leaf share of uptake = C_l VWTLEAF / (VWT + C_l VWTLEAF))

### Photosynthesis Reduction (2)
- `c_leaf_reduc_slope`: Slope of the leaf-salt reduction, per mol/m3 above the threshold (baseline 0.0355, the Campos-Villareal fit). Reduction factor = min(slope * max(c_leaf - threshold, 0), 0.9999)
- `c_leaf_reduc_threshold`: c_leaf where the reduction starts (mol/m3; baseline 12.0)

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
| 1 | `sensitivity_analysis/Phase 1 - pi0/` | `pi0_leaf` and `pi0_stem`, covaried low/med/high | 3 per scenario |
| 2 | `sensitivity_analysis/Phase 2 - root conductance x FE/` | `kr` (baseline and a tenfold reduction) x `E` (0.95, 0.97, 0.98) | 6 per scenario |
| 3 | `sensitivity_analysis/Phase 3 - leaf uptake x photo slope/` | `leaf_uptake_frac` (C_l) x `c_leaf_reduc_slope`, 3 x 3 factorial | 9 per scenario |
| 4+ | to be defined | | |

### Phase 1 - pi0
- Low: `pi0_leaf` = -1.1, `pi0_stem` = -0.9
- Med (baseline): `pi0_leaf` = -1.4, `pi0_stem` = -1.1
- High: `pi0_leaf` = -1.7, `pi0_stem` = -1.3
- Scenario: Drydown (S2) by default; set `active_scenarios` in the config. S2 starts both compartments at `s_init` = 0.42 and dries both after burn-in. Compartment depths are the species rooting depth (`zr_arr: species`, 0.2 m for Pecan), so soil-root conductance and the soil water balance use the same depth. With 0.2 m compartments s falls from 0.42 to about 0.255 over 15 days; transpiration declines from about day 3 and reaches its cuticular floor on days 11-13, so the pi0 runs stay separate for most of the window (starting from 0.32 they merged by day 8; s_init sweep 0.32-0.45). (Earlier 1 m runs: starting from 0.45 barely dried the soil, and starting below about 0.29 put transpiration at its floor from day 0; starting below 0.26 made the leaf water potential solver fail)
- Simulation: 7-day burn-in plus 15-day analysis window (weather is `Lab_Weather_Summer_Average_Day_30min.xlsx`: the time-of-day mean of the lab record, repeated for 45 days)
- Figures: 15-day time series of transpiration, root uptake, soil water potential, storage fluxes, psi_l and flux partition (share of E from roots, stem and leaf storage)
- Final batch: `20261002_1614` (after the gp lag fix in `hydraulics.py`)
- `S2_one` (both sides drying, compartment 1 starting at 0.32) is no longer in the config. S1 (wet, fresh) and S3 (compartment 1 dries, compartment 2 held at 100 mM with salt uptake) stay defined and inactive. A shared scenario list across phases was left for later; each phase still defines its own keys.
- The current config is a 5-day run of S2 and `S2_dry` (7-day burn-in plus 5 days), with `LAI` 8.5 and `VWT` 0.000854. Batch `20261025_2138` still includes the retired `S2_one` runs. `S2_dry` holds compartment 2 at s = 0.32 and dries compartment 1 from s = 0.25. The pi0 runs are group `pi0` on every active scenario. Stem storage at 1x, 5x, and 10x baseline is group `vwt` on `S2_dry` only. On `S2_dry` the three pi0 runs share about 0.03-0.04 mm of nightly efflux into the drying side, because the stem refills above turgor loss. The 5x and 10x volumes cut that efflux to about 0.01-0.02 mm and about 0-0.01 mm. Nightly redistribution figures show compartment 1 only and have no title. Batch `20261025_2105` used s = 0.17 and 0.20 instead. Batch `20261025_1712` is the previous 5-day S2 and `S2_one` result. Batch `20261025_1725` tried `LAI` 2 on that window. The soil only fell from 0.42 to about 0.32, stem relative water content stayed near 0.96, and the pi0 lines stayed together. Nightly efflux into the drier side lasted about four nights. The 14-day test remains batch `20261025_1656`, and the 15-day S2 result remains batch `20261002_1614`.
- Extra figures, written for each active scenario under its run group: relative soil moisture, nightly hydraulic redistribution into compartment 1, the 5-night sum of that efflux as bars, and storage water potential overlaid on soil water potential (solid storage, dashed soil 1, dotted soil 2).
- Live `VWT` is 0.000854. Batch `20261025_1656` reruns S2 and `S2_one` with that volume (`LAI` still 8.5). Both soils still end near s = 0.20 and stem water potential near -4 MPa. The equal-start drydown crosses relative water content 0.90 on day 2.6, as in the `VWT` = 0.0002 batch. The uneven start crosses it on day 2.5 instead of day 1.6.

#### Phase 1 diagnostics
Scripts used to choose the S2 drydown and explain the model behaviour. They are stowed for the record, not part of the standard workflow. Each script's docstring states the conditions it ran under; outputs are written next to the script. Both s_init sweeps and the diagnostics below ran before the gp lag fix, so their E values are somewhat inflated (most for the high run).
- `diagnostics/s_init_sweep/round1_zr1m/`: 1 m compartments. `drydown_sweep.py` (med run, s_init 0.30-0.23) and `drydown_pick.py` (all runs, 0.32/0.30/0.29), which led to s_init 0.32.
- `diagnostics/s_init_sweep/round2_zr02m/`: 0.2 m compartments. `run_sweep.py` (all runs, s_init 0.34-0.45; CSVs gitignored) and `summarize_sweep.py` (table and comparison figure), which led to s_init 0.42.
- `diagnostics/root_conductance/`: soil-root conductance vs soil water, the psi_l collapse at s_init 0.23, and the s = 0.45 vs 0.32 day-1 comparison.
- `diagnostics/night_psi_l/`: why night psi_l stays far below psi_s once the soil dries (stem storage sets it).
- `diagnostics/checks/`: quick checks (compartment depths, batch summary, plot ranges).
- `extra_plots/flux_partition/`: daily supply shares and prototype flux-partition figures (the plot itself is now `flux_partition` in `post_processing.py`).

### Phase 2 - root conductance x FE
- Two run groups (`group` key on each run; figures go to `figures/<scenario>/<group>/`), each a full 3 x 3 factorial:
  - `kr_E`: root radial permeability `kr` = 1e-9, 1e-8 (baseline), 1e-7 m/s x `E` = 0.80, 0.90, 0.95
  - `B_E`: root length density `B_param` = 5000, 10000 (baseline), 20000 x the same `E` levels
- `kr` is now a `HalophyteStemLeafStorageMultiComp` argument (default 1e-8, the value previously hard-coded in `gsr`). `B_param` changes both root surface area and the soil path length; `kr` only the radial term.
- `salt_uptake` and `dynamic_E` on; `E` is the starting efficiency.
- Scenarios (both wet, soil held constant): S4 has 100 mM in compartment 2 only; S5 has 100 mM in both. In S4 the salty compartment's psi_s (about -0.50 MPa) stays below the root node (about -0.29 MPa), so roots lose water to it and take up no salt; E has no effect. S5 forces uptake from salty soil.
- Extra figures: `photosynthesis` (A and the leaf-salt reduction factor, with the reduction onset marked per run) and `leaf_salt` (c_leaf with the 24.3 mol/m3 reduction threshold, and cumulative leaf salt uptake).
- Current batch: `20261004_1928` (S4 + S5), after fixing the initial storage salt mass in `HalophyteStemLeafStorageMultiComp` (it used `LA` instead of `LAI`, so c_stem and c_leaf dropped from 5 to ~0.002 mol/m3 on the first step). In S5, c_leaf rises from 5 to about 7.9 / 6.3 / 5.4 mol/m3 for E = 0.80 / 0.90 / 0.95, so the 24.3 threshold is not reached. Root conductance changes uptake by about 15% (mostly through transpiration). Dynamic E does not engage because c_stem stays far below E x c_max.
- Live baseline `VWT` is 0.0002, the current Pecan value. Batch `20261025_1531` reruns the kr x E factorial on S5 with that volume (the previous S5 batch `20261018_1610` used 0.011). Stem and leaf concentrations rise together and all nine runs reach about 95 and 87 mol/m3 by day 30. Photosynthetic reduction turns on by day 0.5 and daytime photosynthesis is gone by about day 5. `kr` and `E` change the first days, not the 30-day endpoint.
- Batch `20261025_1538` is the same factorial with soil salt at 30 mol/m3 in both compartments. Concentrations level off near the soil value: c_stem about 29-30 mol/m3 and c_leaf about 24-25. The salt reduction factor settles near 0.42-0.44, so daytime photosynthesis stays on (late-window peaks about 6-7 umol/m2/s, down from about 12 on day 1). Onset is still day 0.5-1.4. `kr` and `E` still barely separate the 30-day endpoint.
- Batch `20261025_1558` is a canopy trial on that 30 mol/m3 setup: `LAI` 2 and leaf area 0.0842 m2, with `VWT` raised to 0.000854 so the same stem cylinder is divided by the new leaf area. `species_traits.py` is unchanged. Stem water per ground area stays about 0.0017 m3/m2, so the runs still finish at c_stem about 29-30 and a leaf osmotic shift of about -0.083 to -0.086 MPa. Transpiration per ground area drops about fourfold with `LAI`. That slows the approach to the cap: onset is day 1.5 at `E` 0.85 and day 3.6 at `E` 0.95. `kr` no longer separates the runs, because the lower demand is supplied even at 1e-9 m/s.
- Live `VWT` is 0.000854 and live `LAI` is 8.5, the non-split sapling value. The factorial is `kr` = 1e-9 and 1e-8 (baseline) m/s, `E` = 0.95 (baseline), 0.97, and 0.99, on S5 at 30 mol/m3. That result is batch `20261025_1743`. At baseline kr, onset is day 3.6 / 6.1 / 17.7. The reduced kr adds about a day at 0.95 and 0.97, and reaches day 20.7 at 0.99. The 0.99 runs finish below the soil cap: c_stem about 20-21 and a leaf osmotic shift of about -0.04 MPa, against about 29 and -0.08 MPa at 0.95 and 0.97. Batch `20261025_1657` used 0.98 instead of 0.99. Batch `20261025_1726` tried `LAI` 2 with the 0.98 set.

#### Phase 2 diagnostics
- `diagnostics/salt_mass/salt_mass.py <batch> [scenario]`: salt mass budget per run (water uptake, salt arriving at roots, admitted, stem/leaf split, share of soil salt stock, per-plant mg NaCl) and a cumulative salt mass figure per group. For `20261004_1928` S5: 5.2-6.0 mol/m2 of salt reaches the roots over 15 days, and (1 - E) of it is admitted (implied E equals the set E). That is 0.26-1.2 mol/m2, or 7-33 mg NaCl per plant. Only 0.23% goes to the leaf (`leaf_uptake_frac` 0.5 weighted by VWTLEAF / (VWTLEAF + VWT)); if all admitted salt went to leaf storage, c_leaf would reach 320-1460 mol/m3. The soil's salt mass is held fixed, so uptake (2-8% of the soil stock over 15 days) never depletes it and excluded salt never accumulates.

Candidate parameters for later phases: storage volumes (`VWT`, `VWTLEAF`), storage conductances (`GWMAX`, `GWMAXLEAF`), elastic moduli (`eta_stem`, `eta_leaf`), and salinity settings (`cs_init`, `E`, `dynamic_E`).

### Phase 3 - leaf uptake fraction x photosynthesis reduction slope
- One 3 x 3 factorial. `leaf_uptake_frac` (C_l) = 0.75, 1.0 (baseline), 1.25. `c_leaf_reduc_slope` = 0.026625, 0.0355 (baseline), 0.044375, which is the Campos-Villareal slope plus or minus 25%. The onset threshold stays at 12 mol/m3.
- C_l is the ratio of concentration increases. With `VWT` = 0.011 and `VWTLEAF` = 5e-5 the leaf's share of salt mass is only about 0.34% / 0.45% / 0.57% at those three C_l values, so most admitted salt still goes to the stem.
- The slope enters photosynthesis as min(slope * max(c_leaf - 12, 0), 0.9999). Stomatal conductance is `gsw` = 1.6 * gsc(A) + gcut, and transpiration is computed from `gsw`, so both move only once c_leaf passes 12 and only in proportion to A.
- Scenario: S5, soil held wet. The live config is 30 mol/m3 in both compartments. Baseline otherwise matches Phase 2 (`E` = 0.95, `kr` = 1e-8, salt uptake and dynamic E on). Simulation: 7-day burn-in plus 30-day window. The 100 mol/m3 result is batch `20261025_1659`.
- Figures: transpiration, midday stomatal conductance (daily maximum of gsw), photosynthesis (A and the reduction factor, onset marked), storage salt (c_leaf with the 12 mol/m3 line, and c_stem), pi0 of leaf and stem, leaf water potential, and the salt-response table (30-day leaf delta pi0 and days to reduction onset). Root uptake, soil water potential and storage fluxes are left off: on a wet constant-soil run they copy transpiration or stay flat.
- `c_leaf_reduc_slope` and `c_leaf_reduc_threshold` are constructor arguments on `C3_c_leaf_reduc` (defaults 0.0355 and 12.0). Configs that omit them, including Phase 1, Phase 2 and the lab runs, keep the previous hard-coded fit.
- Batch `20261020_2118` (S5). c_leaf starts at 7 mol/m3 and reaches about 15.8 / 18.5 / 21.0 at C_l = 0.75 / 1.0 / 1.25, so every run crosses 12. Onset is day 17.5 / 13.4 / 10.5 and does not depend on the slope, because the factor is zero until the threshold. By day 30 the reduction factor is 0.10-0.38. Daytime peak transpiration falls from 0.328 um/s before onset to 0.309-0.234 um/s (about 6-29% lower); gsw falls from 0.137 to 0.128-0.095 mol/m2/s and A from 12.3 to 11.2-7.9 umol/m2/s. The widest split is at C_l = 1.25, where c_leaf is furthest above 12, so the slope multiplies a larger excess. Higher slope then slightly slows further salt uptake (less transpiration), so end-of-window c_leaf and leaf delta pi0 are a little smaller on the steep runs. Leaf delta pi0 is about -0.045 / -0.058 / -0.070 MPa. Stem c_stem and pi0_stem barely separate (most salt still goes to the stem; higher C_l makes stem adjustment slightly smaller). Midday psi_l relaxes from about -0.85 MPa toward -0.75 MPa on the strongest run.
- Batch `20261025_1659` reruns that factorial with `VWT` = 0.000854 and soil salt still 100 mol/m3. The smaller stem fills past the threshold immediately: onset is day 0.8-1.5 and the reduction factor is 1 by the end of the window, so daytime photosynthesis is off. End c_leaf is about 40-59 mol/m3. The steeper slope reaches full reduction at a lower concentration, so those runs stop taking up salt sooner and finish with less leaf salt.
- Batch `20261025_1713` is the current 30 mol/m3 factorial, at `LAI` 8.5. Stem concentration still meets near 29 mol/m3 in every run. Leaf concentration separates by C_l into bands near 19, 24, and 28 mol/m3, and the slope barely moves those bands. Onset is day 4.6 / 3.6 / 2.8 at C_l = 0.75 / 1.0 / 1.25. The reduction factor stays below 1, so midday stomatal conductance and photosynthesis stay split: late peaks of gsw are about 0.05-0.12 mol/m2/s and of A about 4-10 umol/m2/s, lowest at high C_l and the steep slope. Batch `20261025_1728` repeated it at `LAI` 2 and matched these concentrations and conductances, because salt mass and storage volume both scale with `LAI`.

---

## 4. Workflow

```
python sensitivity_analysis/sensitivity_analysis.py "sensitivity_analysis/Phase 1 - pi0/config.yaml"
python sensitivity_analysis/post_processing.py "sensitivity_analysis/Phase 1 - pi0/config.yaml" --batch latest
```

- `sensitivity_analysis.py` runs every run x active scenario and saves all outputs. Everything it does is set in the config.
- `post_processing.py` has one function per figure type. It reads the `plots` section of the **current** config, so axis limits and labels can be changed and figures re-made without re-running the model.

### Phase config sections
- `simulation`: species, weather file, timestep, burn-in and analysis days, compartment depths (`zr_arr`: a list in m, or `species` for the species rooting depth in every compartment), `gwmax_mode`. Phase runs keep a uniform depth (0.2 m for Pecan via `zr_arr: species`). Lab experiment configs instead set `ZR_eff = V/A` from measured pot volume and ground area so soil salt concentration matches pot moles / pore volume; see `lab_experiment_simulation/lab_inputs.py`.
- `baseline`: every model parameter, with source comments
- `scenarios` and `active_scenarios`
- `runs`: named overrides of baseline values (add rows for more combinations)
- `plots`: shared settings plus one block per figure

### Folder layout
```
sensitivity_analysis/
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
