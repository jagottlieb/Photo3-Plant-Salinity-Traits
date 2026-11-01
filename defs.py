from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *

class Simulation(object):
	def __init__(self, species_cls, atm_cls, soil_cls, photo_cls, hydro_cls):
		self.species = species_cls
		self.atm = atm_cls
		self.soil = soil_cls
		self.photo = photo_cls
		self.hydro = hydro_cls
	def update(self, dt, phi, ta, qa):
		self.atm.update(phi, ta, qa) 
		self.photo.update(self.atm, self.hydro.psi_l, self.hydro.tl, dt)
		self.hydro.update(self.atm, self.soil, self.photo, dt)
		# Check if soil is SaltySoilGW and pass qgw if needed
		if hasattr(self.soil, 'gw_cs'):  # SaltySoilGW has gw_cs attribute
			self.soil.update(dt, self.species.ZR, self.hydro.qs, self.hydro.qgw)
		else:
			self.soil.update(dt, self.species.ZR, self.hydro.qs)
	def output(self):
		out = {}
		out.update(self.photo.output())
		out.update(self.hydro.output())
		out.update(self.soil.output())
		return out

class SimulationMultiComp(object):
	"""Simulation class for multi-compartment soil with refactored Halophyte hydraulics"""
	def __init__(self, species_cls, atm_cls, soil_cls, photo_cls, hydro_cls, 
	             zr_arr, root_frac_arr, B, dt):
		self.species = species_cls
		self.atm = atm_cls
		self.soil = soil_cls
		self.photo = photo_cls
		self.hydro = hydro_cls
		self.zr_arr = zr_arr
		self.root_frac_arr = root_frac_arr
		self.B = B
		self.dt = dt
		self.num_compartments = len(zr_arr)
		# Initialize storage for qs values (list of lists, one per compartment)
		self.qs_a = [[] for _ in range(self.num_compartments)]
	
	def update(self, dt, phi, ta, qa):
		# Update atmosphere
		self.atm.update(phi, ta, qa)
		
		# Update photosynthesis
		self.photo.update(self.atm, self.hydro.psi_l, self.hydro.tl, dt, soil=self.soil)
		
		# Get salt concentrations if available
		cs_arr = self.soil.cs if hasattr(self.soil, 'cs') and self.soil.cs is not None else None
		
		# Update hydraulics (Halophyte)
		self.hydro.update(
			atm=self.atm,
			soil=self.soil,
			photo=self.photo,
			#s_arr=self.soil.s,
			root_frac_arr=self.root_frac_arr,
			B=self.B,
			cs_arr=cs_arr,
			dt=dt
		)
		
		# Diagnostic: drop hydraulic redistribution so a dry compartment is not rewetted.
		if getattr(self, "clip_negative_qs", False):
			self.hydro.qs = np.maximum(np.atleast_1d(self.hydro.qs), 0.0)

		# Optional: remove salt admitted to the plant from the soil salt stock.
		# Uses the same per-compartment uptake as Halophyte.uptake (positive qs only).
		if (
			getattr(self, "deplete_soil_salt", False)
			and getattr(self.hydro, "Salt_Uptake", False)
			and hasattr(self.soil, "MS")
			and cs_arr is not None
		):
			qs = np.maximum(np.atleast_1d(self.hydro.qs), 0.0)
			cs = np.atleast_1d(cs_arr)
			# Resolve E the same way as Halophyte.uptake (including dynamic_E).
			c_max = self.hydro.c_stem_max if self.hydro.c_stem_max is not None else max(float(np.max(cs)), 1e-6)
			if getattr(self.hydro, "dynamic_E", False):
				c_thresh = self.hydro.E * c_max
				if self.hydro.c_stem < c_thresh:
					E = self.hydro.E
				else:
					E = min(self.hydro.E + (1.0 - self.hydro.E) * (self.hydro.c_stem / c_max), 0.9999)
			else:
				E = self.hydro.E
			# mol/m2 removed this step; match uptake units (um/s * mol/m3 * 1e-6 * dt)
			dMS = qs * cs * (1.0 - E) * 1e-6 * dt
			self.soil.MS = np.maximum(np.atleast_1d(self.soil.MS) - dMS, 0.0)

		# Store qs values for each compartment
		for i in range(self.num_compartments):
			self.qs_a[i].append(self.hydro.qs[i])
		
		# Update soil with new qs from hydraulics
		self.soil.update(dt, self.zr_arr, self.hydro.qs)
	
	def output(self):
		out = {}
		out.update(self.photo.output())
		out.update(self.hydro.output())
		
		# Get soil output - keep compartment data as lists
		soil_out = self.soil.output()
		out.update(soil_out)
		
		# Add qs values
		out['qs'] = self.qs_a
		
		return out

class Atmosphere(object):
	ca = 400. 
	def __init__(self, phi, ta, qa):
		self.phi = phi
		self.ta = ta
		self.qa = qa
		self.cs = self.ca
	def update(self, phi, ta, qa):
		self.phi = phi
		self.ta = ta
		self.qa = qa
