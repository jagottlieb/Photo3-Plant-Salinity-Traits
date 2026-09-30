#!/usr/bin/env python
"""
Systematic Sensitivity Analysis for Photo3-Plant-Salinity-Traits Model

This script performs comprehensive sensitivity analysis using multiple methods:
1. One-at-a-Time (OAT) screening
2. Interaction analysis (2-parameter grids)
3. Sobol' variance-based sensitivity analysis
4. Scenario validation

Author: Sensitivity Analysis Module
Date: 2026-09-29
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Any
import time
import json
from datetime import datetime

# Import model components
from species_traits import *
from soil import SoilMultiple, SaltySoilMultiple, Loam, Sand, Berger, DrydownSoil, ConstantSoil
from defs import SimulationMultiComp, Atmosphere, steps
from photosynthesis import C3, C3_c_leaf_reduc
import hydraulics


class SensitivityAnalyzer:
    """Main class for systematic sensitivity analysis"""
    
    def __init__(self, output_dir: str = "output/sensitivity_analysis"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        (self.output_dir / "oat_results").mkdir(exist_ok=True)
        (self.output_dir / "interaction_heatmaps").mkdir(exist_ok=True)
        (self.output_dir / "sobol_results").mkdir(exist_ok=True)
        (self.output_dir / "parameter_rankings").mkdir(exist_ok=True)
        (self.output_dir / "scenario_results").mkdir(exist_ok=True)
        
        self.results = {}
        self.parameter_space = {}
        
        # Set plot style
        sns.set_style("whitegrid")
        plt.rcParams.update({
            "font.family": "serif",
            "font.serif": ["Caladea", "STIX Two Text", "Constantia", "Georgia", "Cambria"],
            "font.size": 10,
        })
    
    def define_parameter_space(self) -> Dict[str, List[Tuple[str, float]]]:
        """
        Define parameter space with high/medium/low values
        
        Returns:
            Dictionary mapping parameter names to (label, value) tuples
        """
        self.parameter_space = {
            # Stem storage parameters
            "pi0_stem": [("Low", -1.4), ("Baseline", -1.1), ("High", -0.8)],
            "wft_stem": [("Low", 0.8), ("Baseline", 1.0), ("High", 1.2)],
            "wr_stem": [("Low", 0.30), ("Baseline", 0.40), ("High", 0.50)],
            "eta_stem": [("Low", 4.0), ("Baseline", 5.8), ("High", 7.5)],
            "mcap_stem": [("Low", 8.0), ("Baseline", 12.0), ("High", 16.0)],
            "vwi_stem": [("Low", 0.70), ("Baseline", 0.85), ("High", 0.95)],
            "c_stem": [("Low", 3.0), ("Baseline", 5.0), ("High", 7.0)],
            "VWT": [("Low", 0.008), ("Baseline", 0.011), ("High", 0.014)],
            "GWMAX": [("Low", 0.001), ("Baseline", 0.002), ("High", 0.003)],
            
            # Leaf storage parameters
            "pi0_leaf": [("Low", -1.7), ("Baseline", -1.4), ("High", -1.1)],
            "wft_leaf": [("Low", 0.8), ("Baseline", 1.0), ("High", 1.2)],
            "wr_leaf": [("Low", 0.35), ("Baseline", 0.46), ("High", 0.55)],
            "eta_leaf": [("Low", 12.0), ("Baseline", 17.5), ("High", 22.0)],
            "mcap_leaf": [("Low", 8.0), ("Baseline", 12.0), ("High", 16.0)],
            "vwi_leaf": [("Low", 0.75), ("Baseline", 0.90), ("High", 0.98)],
            "VWTLEAF": [("Low", 0.002), ("Baseline", 0.003), ("High", 0.004)],
            "GWMAXLEAF": [("Low", 0.0005), ("Baseline", 0.001), ("High", 0.0015)],
            
            # Soil parameters
            "s_init_0": [("Low", 0.35), ("Baseline", 0.45), ("High", 0.55)],
            "s_init_1": [("Low", 0.35), ("Baseline", 0.45), ("High", 0.55)],
            "cs_init_0": [("Low", 0.0), ("Baseline", 0.0), ("High", 100.0)],
            "cs_init_1": [("Low", 0.0), ("Baseline", 0.0), ("High", 100.0)],
            "root_frac_0": [("Low", 0.3), ("Baseline", 0.5), ("High", 0.7)],
            "root_frac_1": [("Low", 0.3), ("Baseline", 0.5), ("High", 0.7)],
            "B_param": [("Low", 5000), ("Baseline", 10000), ("High", 15000)],
            "F_CAP": [("Low", 0.3), ("Baseline", 0.5), ("High", 0.7)],
            
            # Environmental parameters
            "E": [("Low", 0.85), ("Baseline", 0.95), ("High", 1.0)],
            "dynamic_E": [("Low", False), ("Baseline", False), ("High", True)],
            "salt_uptake": [("Low", False), ("Baseline", False), ("High", True)],
            "leaf_uptake_frac": [("Low", 0.3), ("Baseline", 0.5), ("High", 0.7)],
        }
        
        return self.parameter_space


# Continue with remaining methods...
