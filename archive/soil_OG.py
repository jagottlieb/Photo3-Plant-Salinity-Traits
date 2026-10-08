# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 19:24:34 2021

@author: joshg
"""
from scipy.optimize import *
from dics import *
from functions import *
import numpy as np
import time

#Introduces new classes for use with original Photo3 code.
#SoilGW models a soil with a saturated and unsaturated compartment.
#SaltySoilGW is a subclass of SoilGW which includes an osmotic potential for the unsaturated and groundwater compartments.
#Hydro (need to to rename as HydroGW) models hydraulic systems with an unsaturated compartment
#HalophyteGW models models hydraulic systems with salinity in either or both of the soil compartments.
class LentiscPistacia(object):
	NAME = 'L. pist'
	#PTYPE = C3

	ZR = 2 # Armas, 2008. The depth at which gradient in moisture and salinity begins.
	LAI = 3 #3  #
	GCUT = 0.025 #0.025 #
	GA = 324. # From Jones, 1992 for plant height of 2m and wind speed of 2m/s
	GWMAX = 0.002 #0.002 #.002
	VWT = 0.012 #0.0096*LAI #0.000249 # m3/m2
	GPMAX = 7.5 #
	CAP = 0.15
	VCMAX0 = 37 #57.7
	JMAX0 = 60 #98.5
	PSILA0 = -4.55 #-3. # May Need to consider stomatal conductance redution range
	PSILA1 = -0.3 #-0.5 # May Need to consider stomatal conductance redution range
	RAIW = 10.
	capOn = True

class SoilGW(object):
	EVMAX = 3
	SY = 0.2 #Specific yield for a silt aquifer. This part of the code
	def __init__(self, stype, dynamics, zr, s):
		self.PSI_SS = stype.PSI_SS
		self.B = stype.B
		self.KS = stype.KS
		self.N = stype.N
		self.SH = stype.SH
		self.ZR = zr
		self.s = s
		self.s_a = []
		self.psi_s_a = []
		self.psi_gw_a = []
		self.dynamics = dynamics
		self.rain_amt = 0 # this rain amount is in mm!!
		self.sm_inp = s
	def update(self, dt, zr, qs):
		self.s = self.dynamics.snew(self, dt, zr, qs) #This needs to be updated to accomodate groundwater dynamics.
		#self.s = (dt/(self.N*zr*10.**6)*(-qs - (self.evap(self.s)*1000.)/(24.*60*60)- self.leak(self.s))) + self.s
		self.s_a.append(self.s)
		self.psi_s_a.append(self.psi_s(self.s))
		
	def output(self):
		return {'s': self.s_a, 'psi_s': self.psi_s_a}
	def leak(self, s):
	    """Leakage (um/s) """    
	    try:
	        ans = .11574*self.KS*s**(2.*self.B + 3.)         
	    except OverflowError:
	        ans = 0.
	    return ans                          
	def psi_s(self, s):
	    """Soil Potential (MPa)"""
	    return self.PSI_SS*(s**-self.B)  
	def evap(self, s): 
	    """Soil evaporation rate, per unit ground area (mm/day)"""
	    if s > self.SH:
	        return self.EVMAX*(s - self.SH)/(1. - self.SH)
	    else:
	        return 0.

class SaltySoilGW(SoilGW):
	TS = 293. # soil water temp (K)
	GWTS = 293 # groundwater temp (K)
	IV = 2. # van't hoff coefficient for NaCl
	E = 0.92
	SY = 0.5
	GW_Depth = 3.5 #Depth to subsoil compartment in m
	def __init__(self, stype, dynamics, zr, s, cs, cgw, gw_z): #gw_z = depth to wt, gw_cs = salt concentration in gw, S_y = specific yield
		SoilGW.__init__(self, stype, dynamics, zr, s)
		self.s = s
		self.cs = cs # salt concentration in soil, mol/m3
		self.MS = cs*self.ZR*self.N*s # mass of salt in soil, mol/m2
		self.cs_a = []
		self.gw_z = gw_z #
		self.gw_cs = cgw # Set as constant
		self.S_y = self.SY # Constant
		self.psi_s_a = []
		self.gw_z_a = [] #
		self.gw_cs_a = [] # Constant
		self.psi_gw_osm_a = []
		self.psi_gw_mat_a =[]
		self.psi_s_mat_a = [] 
class DrydownSoil(object):
	def __init__(self):
		pass
	def snew(self, soil, dt, zr, qs):
		return (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s

class ConstantSoil(object):
	def __init__(self):
		pass
	def snew(self, soil, dt, zr, qs):
		return soil.s

class StochasticSoil(object):
	"""takes alpha in cm, lda in 1/d"""
	def __init__(self, alpha, lda):
		self.alpha = alpha
		self.lambda_r = lda
	def rain(self, dt, gamma):
		if np.random.random() > self.lambda_r*dt/(3600.*24.):
			return 0.
		else:
			return np.random.exponential(1./gamma)
	def sLoss(self, soil, dt, zr, qs):
		return (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s
	def snew(self, soil, dt, zr, qs):
		gamma = (soil.N*zr*100.)/self.alpha; #Normalized Depth of Rainfall
		return min(1., self.sLoss(soil, dt, zr, qs) + self.rain(dt, gamma))

class RainSoil(object):
	def __init__(self):
		pass
	def sLoss(self, soil, dt, zr, qs):
		return (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s
	def snew(self, soil, dt, zr, qs):
		"""Takes input of rainfall in mm"""
		return min(1., self.sLoss(soil, dt, zr, qs) + soil.rain_amt/(soil.N*zr*1000.))

class SetSoil(object):
	def __init__(self):
		pass
	def sLoss(self, soil, dt, zr, qs):
		return (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s
	def snew(self, soil, dt, zr, qs):
		return soil.sm_inp

class SaltySoil(object):
	TS = 293. # soil water temp (K)
	IV = 2. # van't hoff coefficient for NaCl
	E = 0.95
	def __init__(self, stype, zr, s, cs):
		self.PSI_SS = stype.PSI_SS
		self.B = stype.B
		self.KS = stype.KS
		self.N = stype.N
		self.SH = stype.SH
		self.ZR = zr
		self.s = s
		self.cs = cs # salt concentration in soil, mol/m3
		self.MS = cs*self.ZR*self.N*s # mass of salt in soil, mol/m2
		self.s_a = []
		self.cs_a = []
		self.rain_amt = 0
		self.sm_inp = s
	def update(self, dt, zr, qs):
		self.s = (dt/(self.N*zr*10.**6)*(-qs - (self.evap(self.s)*1000.)/(24.*60*60)- self.leak(self.s))) + self.s
		self.cs = self.MS/(self.s*self.N*self.ZR) # salt concentration in soil, mol/m3
		self.s_a.append(self.s)
		self.cs_a.append(self.cs)
	def output(self):
		return {'s': self.s_a, 'cs': self.cs_a}
	def psi_s(self, s):
		return self.PSI_SS*(s**-self.B) - self.E*self.cs*R*self.IV*self.TS*10.**(-6.)



class Loam(object):
	PSI_SS = -1.43*10.**-3.
	B = 5.39
	KS = 20.
	N = .45
	SH = .19
	def __init__(self):
		pass

class Sand(object):
	PSI_SS = -.34*10**-3
	B = 4.05
	KS = 200.
	N = .35
	SH = .08
	def __init__(self):
		pass

class SandyLoam(object):
	PSI_SS = -.7*10**-3
	B = 4.9
	KS = 80.
	N = .43
	SH = .14
	def __init__(self):
		pass

class LoamySand(object):
	PSI_SS = -.17*10**-3
	B = 4.38
	KS = 100.
	N = .42
	SH = .08
	def __init__(self):
		pass

class Clay(object):
	PSI_SS = -1.82*10**-3
	B = 11.4
	KS = 1.
	N = .5
	SH = .47
	def __init__(self):
		pass