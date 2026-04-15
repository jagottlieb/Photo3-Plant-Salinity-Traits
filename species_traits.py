from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *
from photosynthesis import *

class Oficu(object):
	"""Opuntia ficus-indica (prickly pear)"""
	NAME = 'O. ficu' # species abbreviation (first letter genus, first four letters species)
	PTYPE = CAM # photosynthetic pathway (C3, C4, or CAM)

	# Plant hydraulic paramters
	ZR = 0.3 # rooting depth (m)
	LAI = 3.5 # leaf area index (-)
	GCUT = 0. # cuticular conductance to water vapor (mm/s)
	#GA = 324.
	GA = 30. # atmospheric conductance per unit ground area (mm/s)
	RAIW = 3. # well-watered root area index (-)
	GPMAX = .4 # maximum plant stem hydraulic conductance (um/MPa/s)

	# Plant water storage parameters (only needed for plant hydraulics with capacitance option)
	GWMAX = .02 # max. conductance between water storage tissue and plant xylem (um/MPa/s)
	VWT = .0113 # maximum water storage depth (m)
	CAP = 0.83 # hydraulic capacitance (MPa-1)

	# Photosynthetic parameters
	VCMAX0 = 18. # maximum carboxylation capacity (umol/m2/sec)
	JMAX0 = 36. # maximum electron transport capacity (umol/m2/sec)
	PSILA0 = -3. # leaf water potential at point of full stomatal closure (MPa)
	PSILA1 = -0.5 # leaf water potential at onset of stomatal closure (MPa)
	
	# CAM-specific parameters (only needed for CAM photosynthetic species)
	MMAX = 230000000.  # max concentration of malic acid (umol/m3) 
	AMMAX = 14. # rate of malic acid storage flux (umol/m2/s) 

class Pmenz(object):
	"""Pseudotsuga menziesii (Douglas fir)"""
	NAME = 'P. menz' # species abbreviation (first letter genus, first four letters species)
	PTYPE = C3 # photosynthetic pathway (C3, C4, or CAM)

	# Plant hydraulic paramters
	ZR = 0.65 # rooting depth (m)
	LAI = 8.4 # leaf area index (-)
	GCUT = .007 # cuticular conductance to water vapor (mm/s)
	GA = 324. # atmospheric conductance per unit ground area (mm/s)
	RAIW = 10. # well-watered root area index (-)
	GPMAX = 0.056 # maximum plant stem hydraulic conductance (um/MPa/s)

	# Plant water storage parameters (only needed for plant hydraulics with capacitance option)
	capOn = True
	GWMAX = .005 # max. conductance between water storage tissue and plant xylem (um/MPa/s)
	VWT = 0.27/LAI # maximum water storage depth (m)
	CAP = 0.15 # hydraulic capacitance (MPa-1)

	# Photosynthetic parameters
	VCMAX0 = 57.7 # maximum carboxylation capacity (umol/m2/sec)
	JMAX0 = 98.5 # maximum electron transport capacity (umol/m2/sec)
	PSILA0 = -3. # leaf water potential at point of full stomatal closure (MPa)
	PSILA1 = -0.5 # leaf water potential at onset of stomatal closure (MPa)
	
	
class Sbico(object):
	"""Sorghum bicolor"""
	NAME = 'S. bico' # species abbreviation (first letter genus, first four letters species)
	PTYPE = C4 # photosynthetic pathway (C3, C4, or CAM)

	# Plant hydraulic paramters
	ZR = 0.5 # rooting depth (m)
	LAI = 5. # leaf area index (-)
	GCUT = 0.1802 # cuticular conductance to water vapor (mm/s)
	GA = 61. # atmospheric conductance per unit ground area (mm/s)
	RAIW = 5.6 # well-watered root area index (-)
	GPMAX = 0.13 # maximum plant stem hydraulic conductance (um/MPa/s)

	# Plant water storage parameters (only needed for plant hydraulics with capacitance option)
	GWMAX = 0. # max. conductance between water storage tissue and plant xylem (um/MPa/s)
	VWT = .000001 # maximum water storage depth (m)
	CAP = 0.15 # hydraulic capacitance (MPa-1)

	# Photosynthetic parameters
	VCMAX0 = 39. # maximum carboxylation capacity (umol/m2/sec)
	JMAX0 = 180. # maximum electron transport capacity (umol/m2/sec)
	PSILA0 = -1.8 # leaf water potential at point of full stomatal closure (MPa)
	PSILA1 = -0.5 # leaf water potential at onset of stomatal closure (MPa)

class Taest(object):
	"""Triticum aestivum (winter wheat)"""
	NAME = 'T. aest' # species abbreviation (first letter genus, first four letters species)
	PTYPE = C3 # photosynthetic pathway (C3, C4, or CAM)

	# Plant hydraulic paramters
	ZR = 0.75 # rooting depth (m)
	LAI = 5. # leaf area index (-)
	GCUT = 0.3 # cuticular conductance to water vapor (mm/s)
	GA = 61. # atmospheric conductance per unit ground area (mm/s)
	RAIW = 5.6 # well-watered root area index (-)
	GPMAX = 11.7  # maximum plant stem hydraulic conductance (um/MPa/s)

	# Plant water storage parameters (only needed for plant hydraulics with capacitance option)
	GWMAX = 0. # max. conductance between water storage tissue and plant xylem (um/MPa/s)
	VWT = .000001 # maximum water storage depth (m)
	CAP = 0.15 # hydraulic capacitance (MPa-1)

	# Photosynthetic parameters
	RD0 = 4.93 # standard dark respiration at 25 C (umol/m2/sec)
	HAV = 62000. # activation energy for Vcmax (J/mol)
	HDV = 202900. # deactivation energy for Vcmax (J/mol)
	VCMAX0 = 83. # maximum carboxylation capacity (umol/m2/sec)
	JMAX0 = 132. # maximum electron transport capacity (umol/m2/sec)
	PSILA0 = -2. # leaf water potential at point of full stomatal closure (MPa)
	PSILA1 = -0.7 # leaf water potential at onset of stomatal closure (MPa)

class Pecan(object): #Pecan tree ()

	NAME = 'Pecan'
	PTYPE = C3

	ZR = 0.25  # rooting depth (m)
	LAI = 1  # leaf area index (m2/m2)
	GCUT = (0.22 * 10**-3 * VW) * 10**3  # Average cuticular rate from Cao, 2019 (mmol/m2/s) converted to mm/s
	GA = 185  # atmospheric conductance (mm/s)
	RAIW = 14.3256  # Woodroof (1934), well-watered root area index (m2/m2)

	# Maximum xylem conductance per unit leaf area (um/MPa/s)
	# Based on Steinberg et al. conversion and scaled by average lab plant measurements.
	GPMAX = 1000 * 8 * 10**-5 * (3900 / ((79 / 4) ** 2) * np.pi) * (((14.6 / 4) ** 2) * np.pi / 1010) / LAI

	# Conductance/storage parameters
	GWMAX = 0.054 * 0.0043 / (0.63 * 4 * 60 * 60)  # max conductance between storage water and xylem (um/MPa/s)
	GWMAXLEAF = 0.0005  # Value from American Beech (um/MPa/s)
	VWT = 0.27 / LAI * 0.1  # max PWS (m3/m2 leaf area); PWS assumed negligible
	VWTLEAF = 0.0504 * (1 / 1000)  # max leaf water storage (m3/m2 leaf area)
	CAP = 0.15  # intrinsic plant hydraulic capacitance (MPa-1)

	RD0 = 3.01  # Standard dark respiration at 25 C (umol/(m^2s))
	HAV = 62000.  # Activation energy for Vc,max (J/mol)
	HDV = 202900.  # Deactivation energy for Vc,max (J/mol)
	VCMAX0 = 80. # maximum carboxylation capacity (umol/m2/sec) Shen et al. Preprint for Pecan Seedlings
	JMAX0 = 175. # maximum electron transport capacity (umol/m2/sec) Shen et al. Preprint for Pecan Seedlings
	# Values taken from Othman et al. 2014
	PSILA0 = -2.1 # leaf water potential at point of full stomatal closure (MPa)
	PSILA1 = -0.8 # leaf water potential at onset of stomatal closure (MPa)

class Pgran(object):
    """Punica granatum (pomegranate)"""
    NAME = 'P. gran' # species abbreviation (first letter genus, first four letters species)
    PTYPE = C3 # photosynthetic pathway (C3, C4, or CAM)

    # Plant hydraulic paramters
    ZR = 0.508 # rooting depth (m) #from USDA PLANTSdata
    LAI = 1.64 # leaf area index (-) #(Meshram et al., 2019)
    GCUT = 0.08 # cuticular conductance to water vapor (mm/s)
    GA = 98.64 # atmospheric conductance per unit ground area (mm/s)
    RAIW = 5.0 # well-watered root area index (-)
    GPMAX = (1/0.01)/LAI # maximum plant stem hydraulic conductance (um/MPa/s)

    # Plant water storage parameters (only needed for plant hydraulics with capacitance option)
    GWMAX = 0.005 # max. conductance between water storage tissue and plant xylem (um/MPa/s)
    VWT = 0.27/LAI # maximum water storage depth (m)
    CAP = 0.4 # hydraulic capacitance (MPa-1)

    # Photosynthetic parameters
    RD0 = 3.01 # standard dark respiration at 25 C (umol/m2/sec)
    HAV = 62000. # activation energy for Vcmax (J/mol)
    HDV = 202900. # deactivation energy for Vcmax (J/mol)
    VCMAX0 = 112.4 # maximum carboxylation capacity (umol/m2/sec) #(Serbin et al., 2015)
    JMAX0 = 224.8 # maximum electron transport capacity (umol/m2/sec) #2 times of VCMAX0
    PSILA0 = -3.47 # leaf water potential at point of full stomatal closure (MPa) #(Rodriguez et al)
    PSILA1 = -1.52 # leaf water potential at onset of stomatal closure (MPa) #(Rodriguez et al)
