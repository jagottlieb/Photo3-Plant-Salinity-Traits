# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 19:24:34 2021

@author: joshg
"""
from scipy.optimize import fsolve
from dics import *
from functions import *
import numpy as np



def F_red(psi_l, c1, c2):  
    """Reduction Function"""
    return exp(-(-psi_l/c1)**c2)

class Atmosphere(object):
	def __init__(self, phi, ta, ca, rh):
		self.phi = phi
		self.ta = ta
		self.ca = ca
		self.qa = qaRh(rh, ta)


class PhotoVico(object):
	KAPPA_2 = .3 # Quantum yield of photosynthesis (mol CO2/mol photon)
	TO = 298 #293.2 # Reference Temperature for photosynthetic parameters (K)
	GAMMA_0 = 34.6 # Parameter for temp dependence of CO2 compensation point (umol/mol) (Leuning, 1995)
	GAMMA_1 = .0451 # Parameter for temp dependence of CO2 compensation point (1/K)
	GAMMA_2 = .000347 # Parameter for temp dependence of CO2 compensation point (1/K^2)
	KC0 = 300 #302. # Michaelis constant for C02 at TO (umol/mol)
	KO0 = 300 #256. # Michaelis constant for 02 at TO (mmol/mol)
	OI = .209  # Oxygen Concentration (mol/mol)
	SVC = 649. # Entropy term for carboxylation (J/mol)
	SVQ = 646. # Entropy term for e-transport (J/mol)
	HKC =  59430. # Activation Energy for Kc (J/mol)
	HKO =  36000. # Activation Energy for Ko (J/mol)
	HKR =  53000. # Activation Energy for Rd (J/mol)
	HDJ = 200000. # Deactivation Energy for Jmax (J/mol)
	HAJ = 50000. # Activation Energy for Jmax (J/mol)
	RD0 = .32 # Standard Dark respiration at 25 C (umol/(m^2s))
	HAV =  72000.  # Activation Energy for Vc,max (J/mol)
	HDV =  200000. # Deactivation Energy for Vc,max (J/mol)
	c1 = 3*10**6 # Reduction function empirical parameter
	c2 = 8 # Reduction function empirical parameter

	def __init__(self, phi, atm):
		self.VCMAX0 = 57
		self.JMAX0 = 98
		self.phi = phi
		self.ared = 1.
		self.light_atten = 1.
		self.tl = atm.ta ####This will need to be updated when tl is made dynamic and solved systematically
	def v_cmax(self, tl, ared):
	    """Maximum carboxylation rate (umol/(m^2s))"""
	    return ared*self.VCMAX0*exp((self.HAV/(R*self.TO))*(1. - self.TO/tl))/(1. + exp((self.SVC*tl - self.HDV)/(R*tl)))
	def k_o(self, tl):
	    """Michaelis-menten coefficient for O2 (mmol/mol)"""
	    return self.KO0*exp(self.HKO/(R*self.TO)*(1. - self.TO/tl))
	def k_c(self, tl):
		"""Michaelis-menten coefficient for CO2 (umol/mol)"""
		return self.KC0*exp(self.HKC/(R*self.TO)*(1. - self.TO/tl))
	def k1(self, phi, tl):###
		"""Photosynthetic Parameter 1 (Vico 2013) (umol/(m^2*s))"""
		return self.jpar(phi,tl)/4
	def k2(self, phi, tl, ared):###
		"""Photosynthetic Parameter 2 (Vico 2013) (umol/mol)"""
		return self.k1(phi,tl)*self.a2(tl, phi)/self.v_cmax(tl, ared)
	def a2(self, tl, phi):###
		"""Photosynthetic Parameter (Vico 2013) (umol/mol)"""
		return self.k_c(tl)*(1 + (self.OI*1000/self.k_o(tl)))
	def gamma(self, tl):
	    """CO2 compensation point (umol/mol)"""
	    return self.GAMMA_0*(1. + self.GAMMA_1*(tl - self.TO) + self.GAMMA_2*(tl - self.TO)**2.);
	def jmax(self, tl):
	    """Max. e- transport rate (umol/(m^2s))"""
	    return self.JMAX0*exp(self.HAJ/(R*self.TO)*(1. - self.TO/tl))/(1. + exp((self.SVQ*tl - self.HDJ)/(R*tl))) 
	def j(self, phi, tl):
	    """Electron transport rate (umol/(m^2s))"""
	    return min((phi*10.**6)/(EP*NA)*self.KAPPA_2*.5, self.jmax(tl)) 
	def jpar(self, phi, tl):
	    """Electron transport rate (umol/(m^2s), based off of PAR, not total solar radiation) Linear for now, upt to Jmax. Vico presents a quadratic formulation."""
	    return min(phi*self.KAPPA_2, self.jmax(tl)) 
	def a_phigsTl(self, gs, rd, k1, k2, ca): ######
	    """Net Photsynthetic Demand for CO2 (Vico,2013) (umol/(m^2s^1))"""
	    return 0.5*(k1 + (k2+ca)*gs - rd - np.sqrt((k1 + (k2-ca)*gs - rd)**2 - 4*gs*(-ca*gs*k2 - k2*rd - k1*self.gamma(self.tl)))) #Equation 1
	def r_d(self, tl):
	    """Dark respiration flux (umol/(m^2s))"""
	    return self.RD0*exp(self.HKR/(R*self.TO)*(1. - self.TO/tl))
	def ciNew(self, ca, gs, rd, k1, k2, psi_l):####
	    """CO2 concentration in mesophyll cytosol (ppm)""" 
	    return ca-(self.an(gs, rd, k1, k2, psi_l)/self.gsc(gs, rd, k1, k2, psi_l, lai, mv, ta, qa))
	def gsc(self, psi_l, lai, mv, ta, qa,gsw, gcut, a): 
	    """Stomatal conductance (Vico 2013 Formulation) to CO2, per unit leaf area (mol/m2/s)"""
	    return max((gsw - gcut)/a, 0) #Equation 2
	def an(self, gs, rd, k1, k2, psi_l, ca, c1, c2): 
		"""Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
		return self.a_phigsTl(gs, rd, k1, k2, ca)*F_red(psi_l, c1, c2)


class SoilGW(object):
	def __init__(self, psi_s):
		self.psi_s = psi_s


class Hydro(object):
	def __init__(self, psi_l, Kmax):
		self.GCUT = 0
		self.psi_l = psi_l
		self.Kmax = Kmax

	def evf(self, psi_l, lai, mv, ta, qa, c1, c2, soil):
	    """Transpiration, per unit ground area (mol/sec/m^2)"""
	    return self.gsrp(psi_l, c1, c2)*(soil.psi_s - self.psi_l)/(mv*lai) #Convert to mol/sec/m^2

	def shf(self, tl, ta, lai):
		"""Sensible heat flux (W/m^2), per unit ground area"""
		return CP_A*RHO_A*self.GA*(tl-ta)/1000.*lai
	def gsrp(self, psi_l, c1, c2):####
	    """Soil-Root-Plant Conductance, per unit ground area (kg/(s-Pa))"""
	    return self.Kmax*exp(-(-psi_l/c1)**c2)
	def gsw(self, psi_l, lai, mv, ta, qa, P_ATM, c1, c2, soil): ###
	    """Stomatal conductance to water, per unit leaf area (mol/m2/sec)"""
	    return self.evf(psi_l, lai, mv, ta, qa, c1, c2, soil)*P_ATM*1000/(VPD(ta, qa)) #Equation 2 


