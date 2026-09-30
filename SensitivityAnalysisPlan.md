# Systematic Sensitivity Analysis Plan

## Overview
A comprehensive sensitivity analysis plan for the Photo3-Plant-Salinity-Traits model to systematically explore how parameter variations affect key outputs.

---

## 1. Parameter Categories

### Stem Storage Parameters (9)
- `pi0_stem`: Stem water potential at full turgor (MPa)
- `wft_stem`: Water fraction at turgor loss point
- `wr_stem`: Water fraction at turgor loss (relative to VWT)
- `eta_stem`: Elastic modulus parameter
- `mcap_stem`: Maximum curvature parameter
- `vwi_stem`: Initial water storage volume fraction
- `c_stem`: Stem storage concentration
- `VWT`: Stem storage volume (m³/m²)
- `GWMAX`: Stem storage conductance (um/MPa/s)

### Leaf Storage Parameters (9)
- `pi0_leaf`: Leaf water potential at full turgor (MPa)
- `wft_leaf`: Water fraction at turgor loss point
- `wr_leaf`: Water fraction at turgor loss (relative to VWTLEAF)
- `eta_leaf`: Elastic modulus parameter
- `mcap_leaf`: Maximum curvature parameter
- `vwi_leaf`: Initial water storage volume fraction
- `c_leaf`: Leaf storage concentration
- `VWTLEAF`: Leaf storage volume (m³/m²)
- `GWMAXLEAF`: Leaf storage conductance (um/MPa/s)

### Soil Parameters (5)
- `s_init_arr`: Initial soil moisture [fraction, fraction]
- `cs_init_arr`: Initial soil salinity [mM, mM]
- `root_frac_arr`: Root fraction per compartment
- `B_param`: Soil moisture potential parameter
- `F_CAP`: Storage fraction available for transpiration

### Environmental Parameters (4)
- `E`: Stomatal conductance scaling factor
- `dynamic_E`: Dynamic stomatal conductance
- `salt_uptake`: Salt uptake enabled
- `leaf_uptake_frac`: Leaf salt uptake fraction

---

## 2. Target Outputs

### Primary Water Status (5)
- `psi_l`: Leaf water potential (MPa)
- `psi_w_stem`: Stem storage water potential (MPa)
- `psi_w_leaf`: Leaf storage water potential (MPa)
- `psi_w_osm_leaf`: Leaf osmotic potential (MPa)
- `psi_w_turgor_leaf`: Leaf turgor pressure (MPa)

### Water Fluxes (4)
- `ev`: Transpiration rate (µmol/m²/s)
- `qs`: Soil-to-plant water flux (µmol/m²/s)
- `qw_stem`: Stem storage-to-xylem flux (µmol/m²/s)
- `qw_leaf`: Leaf storage-to-xylem flux (µmol/m²/s)

### Photosynthesis (2)
- `a`: Photosynthesis rate (µmol/m²/s)
- `gsv`: Stomatal conductance (mol/m²/s)

### Integrated Metrics (9)
- `ev_cum`: Cumulative transpiration (mm)
- `Integrated_transpiration`: Total water loss (µmol)
- `Integrated_soil_flux_comp1`: Root uptake from comp1 (µmol)
- `Integrated_soil_flux_comp2`: Root uptake from comp2 (µmol)
- `Integrated_stem_storage_flux`: Stem contribution (µmol)
- `Integrated_leaf_storage_flux`: Leaf contribution (µmol)
- `Mean_psi_l`: Mean leaf water potential (MPa)
- `Mean_psi_w_stem`: Mean stem water potential (MPa)
- `Mean_psi_w_leaf`: Mean leaf water potential (MPa)

---

## 3. Analysis Methods

### Phase 1: One-at-a-Time (OAT) Screening
- Vary 1 parameter at a time, hold others constant
- Test all stem/leaf storage parameters (18 total)
- Measure sensitivity of key outputs
- **Output**: Parameter ranking
- **Simulations**: ~54

### Phase 2: Interaction Analysis
- Test top 6 parameters from Phase 1
- 2-parameter combinations (3×3 grid each)
- **Output**: Interaction heatmaps
- **Simulations**: ~135

### Phase 3: Sobol' Sensitivity Analysis
- Quantify main effects and interaction effects
- Top 8-10 parameters
- **Output**: First-order (S_i) and total-order (S_Ti) indices
- **Simulations**: ~25

### Phase 4: Validation Scenarios
- Baseline, drydown, drydown+salt, high salinity, high/low storage
- **Output**: Scenario comparisons
- **Simulations**: ~6

---

## 4. Implementation

### Script Structure
```python
# sensitivity_analysis.py
class SensitivityAnalyzer:
    def define_parameter_space(self):  # High/Med/Low values
    def run_oat_analysis(self):       # One-at-a-time
    def run_interaction_analysis(self):  # 2-parameter grids
    def run_sobol_analysis(self):     # Variance-based
    def compute_sensitivity_indices(self):  # S_i, S_Ti
    def generate_report(self):        # Plots and tables
```

### Directory Structure
```
output/sensitivity_analysis/
├── oat_results/
├── interaction_heatmaps/
├── sobol_results/
├── parameter_rankings/
└── scenario_results/
```

---

## 5. Deliverables

1. Parameter ranking by sensitivity to each output
2. Interaction heatmaps
3. Sobol sensitivity indices (S_i, S_Ti)
4. Scenario validation plots
5. Visualization suite for publication
6. Recommendations for calibration/uncertainty quantification

---

## 6. Timeline

| Phase | Duration | Simulations |
|-------|----------|-------------|
| OAT Analysis | 1-2 days | ~54 |
| Interactions | 2-3 days | ~135 |
| Sobol' | 3-4 days | ~25 |
| Validation | 1-2 days | ~6 |
| **Total** | **7-11 days** | **~220** |

---

## 7. Notes

- Each sim: ~5-10 min
- Total runtime: ~20-40 hours
- Use parallel processing for Phase 2 & 3
- Save intermediate results regularly
- Validate extreme parameter combinations

---

*Created: 2026-09-29*