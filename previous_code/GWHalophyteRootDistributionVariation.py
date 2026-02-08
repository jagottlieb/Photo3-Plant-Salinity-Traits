# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 11:02:20 2022

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
	GWMAX = 0.002 #.001   May need to consider this
	GPMAX = 7.5 #
	VWT = 0.012 #0.0096*LAI #0.000249 #
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
	GW_Depth = 5 #Depth to subsoil compartment in m
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
		self.psi_s_osm_a = []
	def update(self, dt, zr, qs, qgw):
		self.gw_z_a.append(self.gw_z)
		self.psi_s_a.append(self.psi_s(self.s))
		self.s_a.append(self.s)
		self.s = (dt/(self.N*zr*10.**6)*(-qs - (self.evap(self.s)*1000.)/(24.*60*60)- self.leak(self.s))) + self.s
		self.cs = self.MS/(self.s*self.N*self.ZR) # salt concentration in soil, mol/m3
		self.gw_z = self.gw_z - dt*60*self.leak(self.s)/10**6 + qgw*60*dt/10**6 ###specificyield not taken into account in the code yet.
		self.gw_cs_a.append(self.gw_cs)
		self.cs_a.append(self.cs)
		self.psi_gw_a.append(self.psi_gw(self.gw_cs))
		self.psi_gw_osm_a.append(self.psi_gw(self.gw_cs))
		self.psi_gw_mat_a.append(self.psi_gw_mat(self.gw_cs))
		self.psi_s_mat_a.append(self.psi_s_mat(self.s))
		self.psi_s_osm_a.append(self.psi_s_osm(self.s))
		
	def output(self):
		return {'s': self.s_a, 'cs': self.cs_a,'gw_z': self.gw_z_a, 'psi_s': self.psi_s_a, 'psi_gw': self.psi_gw_a, 'psi_gw_osm':self.psi_gw_osm_a, 'psi_gw_mat':self.psi_gw_mat_a, 'psi_s_mat': self.psi_s_mat_a, 'psi_s_osm': self.psi_s_osm_a}
	def psi_s(self, s): 
		return self.PSI_SS*(s**-self.B) - self.E*self.cs*R*self.IV*self.TS*10.**(-6.)
	def psi_s_mat(self, s):
		return self.PSI_SS*(s**-self.B)
		"""Soil Matric potential (MPa)"""
	def psi_s_osm(self, s):
		"""Soil Osmostic potential (MPa)"""
		return -self.E*self.cs*R*self.IV*self.TS*10.**(-6.)
	def psi_gw(self, gw_cs):
		"""Groundwater Potential (MPa)""" #Matric potential not significant due to saturation
		return self.PSI_SS - self.E*gw_cs*R*self.IV*self.GWTS*10.**(-6.) - self.GW_Depth*g*RHO_W*10**(-6)
	def psi_gw_mat(self, gw_cs):
		"""Groundwater Matric Potential (MPa)""" #Matric potential not significant due to saturation
		return self.PSI_SS
	def psi_gw_osm(self, gw_cs):
		"""Groundwater Osmotic Potential (MPa)"""
		return -self.E*gw_cs*R*self.IV*self.GWTS*10.**(-6.)



class Hydro(object):
	Plant_h = 3.5 #4 Plant height in meters
	def __init__(self, species):
		self.GPMAX = species.GPMAX
		self.GA = species.GA
		self.gp = 2 #species.GPMAX
		self.GCUT = species.GCUT
		self.RAIW = species.RAIW
		self.zr = species.ZR
		self.lai = species.LAI
		self.psi_l_a = []
		self.gp_a = []
		self.gsv_a = []
		self.tl_a = []
		self.ev_a = []
		self.qgw_a = []
		self.qs_a = []
		self.qbx_a = []
	def rai(self, s):
		"""Root area index (-)"""
		return self.RAIW*s**-self.A_ROOT_l

	def evf(self, photo, phi, ta, psi_l, qa, tl, ci, lai, ared, **kwargs):
	    """Transpiration, per unit ground area (um/sec)"""
	    if self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared, **kwargs) < 0.000000001:
	    	return 0.
	    else:
	    	return max(lai*(1./(self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared, **kwargs)*R*ta/P_ATM*1000000.)+1./(self.GA*1000.))**(-1.)\
	    	*RHO_A/RHO_W*(self.qi(tl, psi_l)-qa), 0.)

	def evfPen(self, photo, phi, ta, psi_l, qa, tl, ci,  lai, ared):
		"""Penman-Monteith transpiration (um/sec)"""
		GAMMA_W = (P_ATM*CP_A)/(.622*LAMBDA_W)
		def delta_s(ta):
			return esat(ta)*(C_SAT*B_SAT)/(C_SAT + ta -273)**2
		def drh(ta, qa):
			return VPD(ta, qa)*.622/P_ATM
		GMGSRATIO = 1.
		gsCAM = self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared)*(1.6*(1.+GMGSRATIO)/(1.6+GMGSRATIO))

		return ((LAMBDA_W*GAMMA_W*self.GA/1000.*RHO_A*drh(ta, qa) + delta_s(ta)*phi)*(R*ta/P_ATM)*gsCAM*1000000.*lai)/ \
		(RHO_W*LAMBDA_W*(GAMMA_W*(self.GA/1000. + (R*ta/P_ATM)*gsCAM*lai) + (R*ta/P_ATM)*gsCAM*lai*delta_s(ta)))
		
	def qi(self, tl, psi_l):
	    """Specific humidity internal to leaf (kg/kg)"""
	    try: 
	        ans =  .622*esat(tl)/P_ATM*exp(psi_l*1000000.*VW/R/tl)
	    except OverflowError:
	        ans = 0.
	    return ans

	def gpf(self, psi_l):  #####Can't have gp = 0
	    """Plant conductance, per unit leaf area (um/(s-MPa))"""
	    if psi_l<-7.5:
	        return 0.00001
	    else:
	        return (-0.52 + 3.45/(1+exp((psi_l+5.32)/(-1.18))))*10**-4*(1/RHO_W)*(1/self.Plant_h)*10**6 #0.05*exp(-(-psi_l/1.5)**2)#2 #self.GPMAX #self.GPMAX*exp(-(-psi_l/1.5)**2) # Curve from Villagrosa 2003 gives LSC in kg/(m*s*MPa), convert to um/(s*MPa)
	def shf(self, tl, ta, lai):
		"""Sensible heat flux (W/m^2), per unit ground area"""
		return CP_A*RHO_A*self.GA*(tl-ta)/1000.*lai
# 	def gsr(self, soil, s, zr): # Using the formulation from Huang et al., 2016
# 	    """Soil-Root Conductance, per unit ground area (um/(s-MPa))"""
# 	    rr = 3*10**-3 #effective root radius (m)
# 	    lam_s = 8000 #4000 #root length density: root length per unit soil volume (m/m^3)
# 	    kr = 10**-9 #10**-9 #root permeability
# 	    B = 2*float(pi)*rr*lam_s #root density: root surface per unit soil volume (m^2/m^3)
# 	    l = 0.53/(float(pi)*B)**0.5
# 	    ks = soil.leak(soil.s)*10**-6/l #convert KS from cm/day to m/s, 1.15741e-7
# 	    return (kr*ks/(kr+ks))*101.9*10**6*B*self.zr #(kr*ks/(kr+ks))*101.9*10**6*B*self.zr #0.15  #convert s^-1 (ie m/s-m) to um/(s-MPa)
# 	def ggwr(self, soil, s, zr): #Using the formulation from Huang et al., 2016
# 	    """Soil-Root Conductance of tap root, per unit ground area (um/(s-MPa))"""
# 	    rr = 3*10**-3 #effective root radius (m)
# 	    lam_gw = 8000 #4000 #root length density: root length per unit soil volume (m/m^3)
# 	    kr = 10**-9 #10**-9 #root permeability
# 	    B = 2*float(pi)*rr*lam_gw #root density: root surface per unit soil volume (m^2/m^3)
# 	    l = 0.53/(float(pi)*B)**0.5 #Length scale characterizing the mean radial distance for the movement of water molecules from bulk soil to root system (Vogel, 2013)
# 	    ks = soil.leak(1)*10**-6/l
# 	    # r = 2 ratio of tap root conductance to lateral root conductance
# 	    return ((kr*ks/(kr+ks))*101.9*10**6)*B*self.zgw #((kr*ks/(kr+ks))*101.9*10**6)*B*self.zr
	def gsw(self, photo, phi, ta, psi_l, qa, tl, ci, ared): 
	    """Stomatal conductance to water, per unit leaf area (mol/m2/sec)"""
	    return photo.gsc(phi, ta, psi_l, qa, tl, ci, ared)*1.6 + (self.GCUT*P_ATM/(1000.*R*ta))


class HalophyteGW(Hydro):
	F_CAP = 0.5
	E = 0.92 # filtration efficiency, unitless
	Salt_Uptake = True
	TS = 293. # soil water temp (K)
	GWTS = 293 # groundwater temp (K)
	IV = 2. # van't hoff coefficient for NaCl


	def __init__(self, species, atm, soil, photo, vwi, cw, GWR, SR):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT #0.000001 
		self.vw = vwi*self.VWT
		self.GWR = GWR
		self.SR = SR
		self.MW = cw*self.vw
		self.CAP = species.CAP
		self.cw = self.MW/self.vw
		self.w = self.vw/self.VWT
		self.zgw = 3.5 # depth to gw table
		#self.psi_l, self.tl= fsolve(self.fBal, (-1, 292.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cm, soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw), soil.gw_cs),  epsfcn = 0.000000001, maxfev=2000) #Played around with xtol and epsfcn
		#self.psi_l, self.tl = self.leaf_solver(atm.ta, self.lai, atm.phi, photo, atm.qa, photo.cx, self.gp, 1, self.vw, self.cw, dt, soil)
		self.out = least_squares(self.fBal, (-1, 292.),args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cm, soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw), soil.gw_cs), method = 'trf', ftol = 3e-16, xtol= 3e-16, x_scale='jac', max_nfev=2000) #least squares solver
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)
		self.qs = self.qsf(soil, soil.s, self.zr, soil.psi_s(soil.s), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw,self.cw), self.gp, self.gwf(self.psi_wf(self.vw,self.cw)), self.lai, self.cw, soil),self.psi_l)
		self.qgw = self.qgwf(soil, soil.psi_gw(soil.gw_cs), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw,self.cw), self.gp, self.gwf(self.psi_wf(self.vw,self.cw)), self.lai, self.cw, soil),soil.gw_cs)
		self.hr_cum = 0
		self.vw_a = []
		self.cw_a = []
		self.psi_b_a = []
		self.psi_x_a = []
		self.qw_a = []
		self.w_a = []
		self.psi_w_a =[]
		self.psi_w_turgor_a = []
		self.psi_w_osm_a =[]
		self.gs_a = []
		self.gwf_a = []
		self.ggwr_a = []
		self.gsr_a = []
		self.hr_cum_a = []
		self.energy_balance = []
		self.water_balance = []
		self.uptake_a =[]
		self.MW_a =[]

	def update(self, atm, soil, photo, dt):
		#[self.psi_l, self.tl], infodict, ier, mesg  = fsolve(self.fBal, (-1, 292.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cm, soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw), soil.gw_cs), epsfcn = 0.0000000000001, full_output=True, maxfev = 2000) #Played around with xtol and epsfcn. Output infodict to diagnose issues
		#print(infodict)
		#print(ier)
		#print(mesg)
		#self.psi_l, self.tl = self.leaf_solver(atm.ta, self.lai, atm.phi, photo, atm.qa, photo.cx, self.gp, 1, self.vw, self.cw, dt, soil)
		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)
		self.MW = self.MW + self.uptake(self.E, self.gp,self.psi_x(self.ev, self.psi_l, self.gp, self.lai), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw, self.cw), self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai, self.cw, soil), self.lai, soil.gw_cs)
		self.uptake_a.append(self.uptake(self.E, self.gp,self.psi_x(self.ev, self.psi_l, self.gp, self.lai), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw, self.cw), self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai, self.cw, soil), self.lai, soil.gw_cs))
		self.cw = (self.MW/self.vw) 
		self.w = self.vw/self.VWT
		self.gp = self.gpf(self.psi_l)
		self.out = least_squares(self.fBal, (-1, 292.),args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cm, soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw), soil.gw_cs), method ='trf', ftol = 3e-16, xtol = 3e-16, x_scale ='jac', max_nfev=2000)
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		#print(self.out.fun)
		#print(self.tl)
		self.energy_balance.append(atm.phi - self.shf(self.tl, atm.ta, self.lai) - LAMBDA_W*RHO_W*self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1)/1000000.)
		self.water_balance.append(self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.) - self.qbx(self.gp,self.psi_x(self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.), self.psi_l, self.gp,self.lai), self.psi_b(self.vw,self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.), self.psi_l, self.psi_wf(self.vw,self.cw),self.gp,self.gwf(self.psi_wf(self.vw,self.cw)), self.lai,self.cw, soil), self.lai) - self.qwf(self.vw,self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.), self.gp,self.psi_l, self.lai, self.cw,dt))
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)
		self.qs = self.qsf(soil, soil.s, self.zr, soil.psi_s(soil.s), self.psi_b(self.vw, self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.), self.psi_l, self.psi_wf(self.vw,self.cw), self.gp, self.gwf(self.psi_wf(self.vw,self.cw)), self.lai, self.cw, soil), self.psi_l)
		self.qgw = self.qgwf(soil, soil.psi_gw(soil.gw_cs), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw,self.cw), self.gp, self.gwf(self.psi_wf(self.vw,self.cw)), self.lai, self.cw, soil),soil.gw_cs)
		self.hr_cum = (self.hr_cum + self.hr(self.qs)*30*60/1000) #Units of mm, updating cumulatively at each 30 min timestep
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, 1.))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
		self.vw_a.append(self.vw)
		self.cw_a.append(self.cw)
		self.qgw_a.append(self.qgw)
		self.gs_a.append(self.gsr(soil, soil.s, self.zr, self.psi_l))
		self.qs_a.append(self.qs)
		self.psi_b_a.append(self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw, self.cw), self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai, self.cw, soil))
		self.psi_x_a.append(self.psi_x(self.ev, self.psi_l, self.gp, self.lai))
		self.qw_a.append(self.qwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt))
		self.qbx_a.append(self.qbx(self.gp, self.psi_x(self.ev, self.psi_l,self.gp, self.lai), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw, self.cw),self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai,self.cw, soil), self.lai))
		self.psi_w_a.append(self.psi_wf(self.vw, self.cw))
		self.psi_w_turgor_a.append(self.psi_wf_turgor(self.vw, self.cw))
		self.psi_w_osm_a.append(self.psi_wf_osm(self.vw, self.cw))
		self.gwf_a.append(self.gwf(self.psi_wf(self.vw, self.cw)))
		self.ggwr_a.append(self.ggwr(soil, soil.s, self.zr, self.psi_l))
		self.gsr_a.append(self.gsr(soil, soil.s, self.zr, self.psi_l))
		self.hr_cum_a.append(self.hr_cum)
		self.w_a.append(self.w)
		self.MW_a.append(self.MW)
		#self.uptake_a.append((self.uptake(self.E, self.gp,self.psi_x(self.ev, self.psi_l, self.gp, self.lai), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw, self.cw), self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai, self.cw, soil), self.lai, soil.gw_cs)/self.vw))
		


	def output(self):
		return {'psi_l': self.psi_l_a, 'psi_w': self.psi_w_a, 'psi_w_osm': self.psi_w_osm_a, 'psi_w_turgor': self.psi_w_turgor_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 'ev': self.ev_a, 'ev_cum': np.cumsum(list(i*1.8 for i in self.ev_a)), 'vw': self.vw_a, 'cw': self.cw_a, 'qgw': self.qgw_a, 'qs': self.qs_a, 'psi_b': self.psi_b_a, 'psi_x': self.psi_x_a, 'qw': self.qw_a, 'qbx': self.qbx_a, 'gs': self.gs_a, 'gwf': self.gwf_a, 'ggwr': self.ggwr_a, 'gsr': self.gsr_a, 'Energy balance': self.energy_balance, 'Water balance': self.water_balance, 'hr_cum': self.hr_cum_a, 'w':self.w_a, 'Uptake':self.uptake_a, 'MW':self.MW_a}

	def psi_wf(self, vw, cw): 
	    """Water potential of stored water (MPa)"""
	    TL = 293.
	    aF = 0.75 #Apoplastic fraction, taken from DeCaceres, 2021
	    psi_wf_osm_ft = -(self.MW/self.VWT)*R*self.IV*TL*10.**(-6.) #Osmotic Potential at full turgor
	    eta = 27.7 #Modulus of Elasticity for Lentiscus following Christofferson, 2016
	    #print(psi_wf_osm_ft)
	    return psi_wf_osm_ft/(vw/self.VWT) + max(-psi_wf_osm_ft - eta*((1-vw/self.VWT)/(1-aF)),0) + self.F_CAP*self.Plant_h*g*RHO_W*10**(-6) #(vw/self.VWT-self.D1)**self.D2 - cw*R*self.IV*TL*10.**(-6.) + 0.5*self.Plant_h*g*RHO_W*10**(-6)#Turgor pressure minus osmotic potential in the storage plus gravitational potential at 0.5 plant height
	def psi_wf_turgor(self, vw, cw): 
	    """Turgor pressure of stored water (MPa)"""
	    TL = 293.
	    eta = 27.7
	    aF = 0.75
	    psi_wf_osm_ft = -(self.MW/self.VWT)*R*self.IV*TL*10.**(-6.)
	    return max(-psi_wf_osm_ft - eta*(1-(vw/self.VWT)),0) #(vw/self.VWT-self.D1)**self.D2
	def psi_wf_osm(self, vw, cw):
	    """Osmotic potential of stored water (MPa)"""
	    TL = 293.
	    aF = 0.75
	    psi_wf_osm_ft = -(self.MW/self.VWT)*R*self.IV*TL*10.**(-6.)
	    return psi_wf_osm_ft/(vw/self.VWT) #-cw*R*self.IV*TL*10.**(-6.)
	def a(self, lai, psi_s, psi_gw, gw, psi_w, soil, s, zr, psi_l): #Algebraic simplification term for hydraulic equations
		return (self.gsr(soil, s, zr, psi_l)*psi_s + self.ggwr(soil, s, zr, psi_l)*psi_gw)
	def b(self, gp, gw, lai): #Algebraic simplification term for hydraulic equations
		return ((self.F_CAP/(gp*lai)) + ((1-self.F_CAP)/(gp*lai)) + ((self.F_CAP*(1-self.F_CAP)*gw/((gp**2)*lai))))
	def c(self, gp, gw, lai): #Algebraic simplification term for hydraulic equations
		return (self.F_CAP*gw/gp)
	def d(self, soil, s, zr, psi_l): #Algebraic simplification term for hydraulic equations
		return (self.gsr(soil, s, zr, psi_l) + self.ggwr(soil, s, zr, psi_l))
	def e(self, gp, lai): #Algebraic simplification term for hydraulic equations
		return (gp*lai/self.F_CAP)
	def psi_x(self, ev, psi_l, gp, lai): 
	    return (ev*(1-self.F_CAP)/(lai*gp) + psi_l)
	def psi_b(self, vw, ev, psi_l, psi_w, gp, gw, lai, cw, soil):
	    return (self.e(gp,lai)*self.psi_x(ev, psi_l, gp, lai) + self.a(self.lai, soil.psi_s(soil.s), soil.psi_gw(soil.gw_cs), self.gwf(self.psi_wf(self.vw, self.cw)), self.psi_wf(self.vw, self.cw), soil, soil.s, self.zr, psi_l))/(self.e(gp,lai)+self.d(soil, soil.s,self.zr, psi_l)) #Update psi_b - it shouldn't be lower than psi_l at peak transpiration...
	def qwf(self, vw, ev, gp, psi_l, lai, cw, dt):
	    """Stored water flux, per unit ground area (um/s)""" 
	    return (vw - self.vwf(vw, ev, gp, psi_l, lai, cw, dt))*lai*10.**6/dt
	def qsf(self, soil, s, zr, psi_s, psi_b, psi_l):
	    """Soil water flux, per unit ground area (um/s)"""
	    return self.gsr(soil, s, zr, psi_l)*(psi_s-psi_b)
	def hr(self, qsf):
		if qsf < 0:
			return -qsf
		else:
			return 0
	def qgwf(self, soil, psi_gw, psi_b, gw_cs):
	    return self.ggwr(soil, soil.s, self.zr, self.psi_l)*(psi_gw - psi_b)
	def qbx(self, gp, psi_x, psi_b, lai):
	    return (gp*lai/self.F_CAP)*(psi_b - psi_x)
	def gwf(self, psi_w):
	    """Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
	    return self.GWMAX*(self.vw/self.VWT)**4 #self.GWMAX*exp(-(-psi_w/2.)**2.) self.GWMAX
	def gsr(self, soil, s, zr, psi_l): # Using the formulation from Huang et al., 2016
	    """Soil-Root Conductance, per unit ground area (um/(s-MPa))"""
	    rr = 0.2*10**-3 #3*10**-3 #effective root radius (m)
	    B = self.SR #Soil default fraction of 0.88 1000 used in Vogel paper #root length density: root length per unit soil volume (m/m^3)
	    kr = 10**-8 #10**-9 #root permeability
	    Ar = 2*float(pi)*rr*B #root density: root surface per unit soil volume (m^2/m^3)
	    l = 0.53/(float(pi)*B)**0.5 #m
	    ks = soil.leak(soil.s)*10**-6/l #to s^-1
	    return (kr*ks/(kr+ks))*101.9*10**6*Ar*self.zr #exp(-(-psi_l/1.5)**1)*(kr*ks/(kr+ks))*101.9*10**6*B*self.zr #(kr*ks/(kr+ks))*101.9*10**6*B*self.zr #(kr*ks/(kr+ks))*101.9*10**6*B*self.zr #0.15  #convert s^-1 (ie m/s-m) to um/(s-MPa)
	def ggwr(self, soil, s, zr, psi_l): #Using the formulation from Huang et al., 2016
	    """Soil-Root Conductance of tap root, per unit ground area (um/(s-MPa))"""
	    rr = 0.2*10**-3 #0.005 cm in Vogel #3*10**-3 #effective root radius (m)
	    B = self.GWR # GW Fraction default of 0.12 4000 #root length density: root length per unit soil volume (m/m^3)
	    kr = 10**-8 #10**-9 #root permeability
	    Ar = 2*float(pi)*rr*B #root surface density: root surface per unit soil volume (m^2/m^3)
	    l = 0.53/(float(pi)*B)**0.5 #Length scale characterizing the mean radial distance for the movement of water molecules from bulk soil to root system (Vogel, 2013)
	    ks = soil.leak(1)*10**-6/l
	    # r = 2 ratio of tap root conductance to lateral root conductance
	    return ((kr*ks/(kr+ks))*101.9*10**6)*Ar*(self.zgw-self.zr) #exp(-(-psi_l/1.5)**1)*((kr*ks/(kr+ks))*101.9*10**6)*B*self.zgw #((kr*ks/(kr+ks))*101.9*10**6)*B*self.zgw 
	def vwf(self, vw, ev, gp, psi_l, lai, cw, dt): 
	    """Stored water volume, per unit leaf area (m3/m2)"""
	    psi_w = self.psi_wf(vw, cw) # Problems at low storage
	    return min(vw - self.gwf(psi_w)*(psi_w - (ev*(1. - self.F_CAP))/(lai*gp) - psi_l)*dt/10.**6, self.VWT) #max(min(vw - self.gwf(psi_w)*(psi_w - (ev*(1. - self.F_CAP))/(lai*gp) - psi_l)*dt/10.**6, self.VWT),0.0001)
	def uptake(self, E, gp, psi_x, psi_b, lai, c_gw):
	    if self.Salt_Uptake == True:
		    return self.qbx(gp, psi_x, psi_b, lai)*c_gw*(1-E)*30*60*10**(-6)
		    #return self.qwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)*c_gw*(1-E)*30*60*10**(-6)
		    print(self.qbx(gp, psi_x, psi_b, lai)*c_gw*(1-E)*30*60)
	    else:
		    return 0
		    print("0")
	def fBal(self, params, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, psi_w, gw_cs): 
	    psi_l, tl = params
	    if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
		    #return (phi*lai - self.shf(tl, ta, lai) - LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)-self.qbx(gp,self.psi_x(self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,gp,lai),self.psi_b(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,psi_w,gp,self.gwf(psi_w),lai,self.cw,soil),lai) - self.qwf(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),gp,psi_l,lai, self.cw,dt))\
		    return (self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)-self.qbx(gp,self.psi_x(self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,gp,lai),self.psi_b(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,psi_w,gp,self.gwf(psi_w),lai,self.cw,soil),lai) - self.qwf(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),gp,psi_l,lai, self.cw,dt), phi*lai - self.shf(tl, ta, lai) - LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000.)
			    #self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
    			#-((self.a(self.lai, soil.psi_s(s), soil.psi_gw(soil.gw_cs), self.gwf(self.psi_wf(self.vw, self.cw)), self.psi_wf(self.vw, self.cw), soil, soil.s, self.zr) + \
    			#(self.d(soil, soil.s, self.zr) + self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)*self.psi_wf(self.vw, self.cw) - \
				#(self.d(soil, soil.s, self.zr)*(1 + self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)) + self.gwf(self.psi_wf(self.vw, self.cw))*self.lai)*psi_l)/ \
				#(1 + (self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)*((1/self.F_CAP) - 1)) + self.d(soil, soil.s, self.zr)*self.b(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)))))
	    else:
	    	return (phi - self.shf(tl, ta, lai) - LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) - self.qbx(gp,self.psi_x(self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,gp,lai), self.psi_b(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,psi_w,gp,self.gwf(psi_w),lai,self.cw,soil),lai) - self.qwf(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),gp,psi_l,lai, self.cw,dt))
		    	#self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) \
				#-((self.a(self.lai, soil.psi_s(s), soil.psi_gw(soil.gw_cs), self.gwf(self.psi_wf(self.vw, self.cw)), self.psi_wf(self.vw, self.cw), soil, soil.s, self.zr) + \
    			#(self.d(soil, soil.s, self.zr) + self.e(self.gp, self.lai))*self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)*self.psi_wf(self.vw, self.cw) - \
				#(((self.d(soil, soil.s, self.zr) + self.e(self.gp, self.lai))*(1 + (self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai))) - self.e(self.gp, self.lai))*psi_l))/ \
				#((-(1-self.F_CAP)/self.F_CAP) + (self.d(soil, soil.s, self.zr) + self.e(sel