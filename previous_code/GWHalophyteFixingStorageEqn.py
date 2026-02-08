# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 19:24:34 2021

@author: joshg
"""
from scipy.optimize import * #fsolve
from dics import *
from functions import *
import numpy as np

class SoilGW(object):
	EVMAX = 3
	SY = 0.2 ####Specific yield for a silt aquifer.
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
		self.s = self.dynamics.snew(self, dt, zr, qs)
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
	E = 0.95
	SY = 0.5
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
		self.psi_s_mat_a = [] 
		self.psi_s_osm_a = []
	def update(self, dt, zr, qs, qgw):
		self.gw_z_a.append(self.gw_z)
		self.psi_s_a.append(self.psi_s(self.s))
		self.s_a.append(self.s)
		self.s = (dt/(self.N*zr*10.**6)*(-qs - (self.evap(self.s)*1000.)/(24.*60*60)- self.leak(self.s))) + self.s
		self.cs = self.MS/(self.s*self.N*self.ZR) # salt concentration in soil, mol/m3
		self.gw_z = self.gw_z - dt*60*self.leak(self.s)/10**6 + qgw*60*dt/10**6 ### specific yield not taken into account 
		self.gw_cs_a.append(self.gw_cs) ##########
		self.cs_a.append(self.cs)
		self.psi_gw_a.append(self.psi_gw(self.gw_cs))
		self.psi_s_mat_a.append(self.psi_s_mat(self.s))
		self.psi_s_osm_a.append(self.psi_s_osm(self.s))
		
	def output(self):
		return {'s': self.s_a, 'cs': self.cs_a,'gw_z': self.gw_z_a, 'psi_s': self.psi_s_a, 'psi_gw': self.psi_gw_a, 'psi_s_mat': self.psi_s_mat_a, 'psi_s_osm': self.psi_s_osm_a}
	def psi_s(self, s): 
		return self.PSI_SS*(s**-self.B) - self.E*self.cs*R*self.IV*self.TS*10.**(-6.)
	def psi_s_mat(self, s):
		return self.PSI_SS*(s**-self.B) 
	def psi_s_osm(self, s): ################
		return -self.E*self.cs*R*self.IV*self.TS*10.**(-6.)
	def psi_gw(self, gw_cs): 
		return -self.E*gw_cs*R*self.IV*self.GWTS*10.**(-6.)



class Hydro(object):
	A_ROOT_l = 8.
	def __init__(self, species, lam_s1, lam_gw1):
		self.GPMAX = species.GPMAX
		self.GA = species.GA
		self.gp = species.GPMAX
		self.GCUT = 0 #species.GCUT
		self.RAIW = species.RAIW
		self.zr = species.ZR
		self.lam_s = lam_s1
		self.lam_gw = lam_gw1
		self.lai = 8.4 # 1 species.LAI
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
	    if self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared, **kwargs) < 0.00001:
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

	    #return .622*esat(tl)/P_ATM*exp(psi_l*1000000.*VW/R/tl)
	def gpf(self, psi_l):  #####Can't have gp = 0
	    #Plant conductance, per unit leaf area (um/(s-MPa))
	    if psi_l<-10:
	        return 0.00001
	    else:
	        return 12*exp(-(-psi_l/1.5)**2) #self.GPMAX*exp(-(-psi_l/2.)**2.) #6*exp(-(-psi_l/2)**1.5) 
	        #return self.GPMAX
	def shf(self, tl, ta, lai):
		"""Sensible heat flux (W/m^2), per unit ground area"""
		return CP_A*RHO_A*self.GA*(tl-ta)/1000.*lai
	def gsr(self, soil, s, zr, lam_s):
	    """Soil-Root Conductance, per unit ground area (um/(s-MPa))"""
	    """return (soil.leak(s)*sqrt(self.rai(s))*1000000.)/(float(pi)*g*RHO_W*zr)"""
	    rr = 3*10**-3 #effective root radius (m)
	    #lam = 5000 #500 1500 4000 root length density: root length per unit soil volume (m/m^3)
	    kr = 10**-9 #root permeability
	    B = 2*float(pi)*rr*lam_s #root density: root surface per unit soil volume (m^2/m^3)
	    l = 0.53/(float(pi)*B)**0.5
	    ks = soil.leak(soil.s)*10**-6/l #convert KS from cm/day to m/s, 1.15741e-7
	    #return self.GPMAX Constant for now
	    return (kr*ks/(kr+ks))*101.9*10**6 #convert s^-1 to um/(s-MPa)
	def ggwr(self, soil, s, zr, lam_gw): ########## 
	    """Soil-Root Conductance of tap root, per unit ground area (um/(s-MPa))"""
	    rr = 3*10**-3 #effective root radius (m)
	    lam = 1500 #500 4000 5000 #root length density: root length per unit soil volume (m/m^3)
	    kr = 10**-9 #root permeability
	    B = 2*float(pi)*rr*lam_gw #root density: root surface per unit soil volume (m^2/m^3)
	    l = 0.53/(float(pi)*B)**0.5
	    ks = soil.leak(1)*10**-6/l
	    # r = 2 ratio of tap root conductance to lateral root conductance
	    return (kr*ks/(kr+ks))*101.9*10**6
	def gsw(self, photo, phi, ta, psi_l, qa, tl, ci, ared): 
	    """Stomatal conductance to water, per unit leaf area (mol/m2/sec)"""
	    #return gsN(phi, Ta, psi_l, qa, Tl, ci, t)*(1.6*(1. + GMGSRATIO[pType[species]]))/(1.6 + GMGSRATIO[pType[species]]) + (gcut[species]*po/(1000.*R*Ta))
	    return photo.gsc(phi, ta, psi_l, qa, tl, ci, ared)*1.6 + (self.GCUT*P_ATM/(1000.*R*ta))


class HalophyteGW(Hydro):
	F_CAP = 0.5
	CW = 150. # default salt concentration in plant, mmol/L
	E = 0.95 # filtration efficiency, unitless
	TS = 293. # soil water temp (K)
	GWTS = 293 # groundwater temp (K)
	IV = 2. # van't hoff coefficient for NaCl
	#CS = 150. # salt concentration in soil, mol/m3
	D1 = .028 # parameter for turgor pressure, for agave
	D2 = 8. # parameter for turgor pressure, for agave

	def __init__(self, species, atm, soil, photo, vwi, cw, lam_s1, lam_gw1):
		Hydro.__init__(self, species, lam_s1, lam_gw1)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.vw = vwi*self.VWT
		self.MW = cw*self.vw
		self.CAP = species.CAP
		self.cw = self.MW/self.vw
		self.psi_l, self.tl = fsolve(self.fBal, (-2, 305.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cm, soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw), soil.gw_cs),xtol=1e-20)
		self.gp = self.gpf(self.psi_l)
		self.lam_s = lam_s1
		self.lam_gw = lam_gw1
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, 1.)
		self.qs = self.qsf(soil, soil.s, self.zr, soil.psi_s(soil.s), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw,self.cw), self.gp, self.gwf(self.psi_wf(self.vw,self.cw)), self.lai, self.cw, soil))
		self.qgw = self.qgwf(soil, soil.psi_gw(soil.gw_cs), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw,self.cw), self.gp, self.gwf(self.psi_wf(self.vw,self.cw)), self.lai, self.cw, soil),soil.gw_cs)
		self.vw_a = []
		self.cw_a = []
		self.psi_b_a = []
		self.psi_x_a = []
		self.qw_a = []
		self.psi_w_a =[]
		self.gs_a = []
		self.gwf_a = []
		self.ggwr_a = []
		self.gsr_a = []
		self.energy_balance = []
		self.water_balance = []

	def update(self, atm, soil, photo, dt):
		self.psi_l, self.tl = fsolve(self.fBal, (-2, 305.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cm, soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw), soil.gw_cs),xtol=1e-20)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, 1.)
		self.gp = self.gpf(self.psi_l)
		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)
		self.cw = self.MW/self.vw
		self.qs = self.qsf(soil, soil.s, self.zr, soil.psi_s(soil.s), self.psi_b(self.vw, self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, 1.), self.psi_l, self.psi_wf(self.vw,self.cw), self.gp, self.gwf(self.psi_wf(self.vw,self.cw)), self.lai, self.cw, soil))
		self.qgw = self.qgwf(soil, soil.psi_gw(soil.gw_cs), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw,self.cw), self.gp, self.gwf(self.psi_wf(self.vw,self.cw)), self.lai, self.cw, soil),soil.gw_cs)
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, 1.))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
		self.vw_a.append(self.vw)
		self.cw_a.append(self.cw)
		self.qgw_a.append(self.qgw)
		self.gs_a.append(self.gsr(soil, soil.s, self.zr, self.lam_s))
		self.qs_a.append(self.qs)
		self.psi_b_a.append(self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw, self.cw), self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai, self.cw, soil))
		self.psi_x_a.append(self.psi_x(self.ev, self.psi_l, self.gp, self.lai))
		self.qw_a.append(self.qwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt))
		self.qbx_a.append(self.qbx(self.gp, self.psi_x(self.ev, self.psi_l,self.gp, self.lai), self.psi_b(self.vw, self.ev, self.psi_l, self.psi_wf(self.vw, self.cw),self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai,self.cw, soil), self.lai))
		self.psi_w_a.append(self.psi_wf(self.vw, self.cw))
		self.gwf_a.append(self.gwf(self.psi_wf(self.vw, self.cw)))
		self.ggwr_a.append(self.ggwr(soil, soil.s, self.zr, self.lam_gw))
		self.gsr_a.append(self.gsr(soil, soil.s, self.zr, self.lam_s))
		self.energy_balance.append(atm.phi*self.lai - self.shf(self.tl, atm.ta, self.lai) -LAMBDA_W*RHO_W*self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)/1000000.)
		self.water_balance.append(self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.) - self.qbx(self.gp,self.psi_x(self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.),self.psi_l,self.gp,self.lai), self.psi_b(self.vw,self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.),self.psi_l,self.psi_wf(self.vw,self.cw),self.gp,self.gwf(self.psi_wf(self.vw,self.cw)),self.lai,self.cw,soil),self.lai) - self.qwf(self.vw,self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.),self.gp,self.psi_l,self.lai,self.cw,dt))
		

	def output(self):
		return {'psi_l': self.psi_l_a,'psi_w': self.psi_w_a,'gp': self.gp_a,'gsv': self.gsv_a,'tl': self.tl_a,'ev': self.ev_a,'vw': self.vw_a,'cw': self.cw_a,'qgw': self.qgw_a,'qs': self.qs_a,'psi_b':self.psi_b_a,'psi_x':self.psi_x_a,'qw':self.qw_a,'qbx':self.qbx_a,'gs':self.gs_a,'gwf':self.gwf_a,'ggwr':self.ggwr_a,'gsr':self.gsr_a,'Energy balance':self.energy_balance,'Water balance':self.water_balance}

	def psi_wf(self, vw, cw): 
	    """Water potential of stored water (MPa)"""
	    TL = 293.
	    return (vw/self.VWT-self.D1)**self.D2 - cw*R*self.IV*TL*10.**(-6.)
	def a(self, lai, psi_s, psi_gw, gw, psi_w, soil, s, zr, lam_s, lam_gw): ###Intake variables for functions
		return (self.gsr(soil, s, zr, lam_s)*psi_s + self.ggwr(soil, s, zr, lam_gw)*psi_gw)
	def b(self, gp, gw, lai): ###Intake variables for functions
		return ((self.F_CAP/(gp*lai)) + ((1-self.F_CAP)/(gp*lai)) + ((self.F_CAP*(1-self.F_CAP)*gw/((gp**2)*lai))))
	def c(self, gp, gw, lai): ###Added
		return (self.F_CAP*gw/gp)
	def d(self, soil, s, zr): ###Intake variables for functions
		return (self.gsr(soil, s, zr, self.lam_s) + self.ggwr(soil, s, zr, self.lam_gw))
	def e(self, gp, lai):
		return (gp*lai/self.F_CAP)
	def psi_x(self, ev, psi_l, gp, lai): 
	    return (ev*(1-self.F_CAP)/(lai*gp) + psi_l)
	def psi_b(self, vw, ev, psi_l, psi_w, gp, gw, lai, cw, soil): ###### Fixed
	    #return (ev + gw*lai*(self.psi_x(ev, psi_l, gp, lai)-psi_w) - self.gsr(soil, soil.s, self.zr, self.lam_s)*soil.psi_s(soil.s)- self.ggwr(soil, soil.s, self.zr, self.lam_gw)*soil.psi_gw(soil.gw_cs))/-self.d(soil, soil.s, self.zr)
	    return (self.e(gp,lai)*self.psi_x(ev, psi_l, gp, lai) + self.a(self.lai, soil.psi_s(soil.s), soil.psi_gw(soil.gw_cs), self.gwf(self.psi_wf(self.vw, self.cw)), self.psi_wf(self.vw, self.cw), soil, soil.s, self.zr, self.lam_s, self.lam_gw))/(self.e(gp,lai)+self.d(soil, soil.s,self.zr)) #Update psi_b - it shouldn't be lower than psi_l at peak transpiration...
	def qwf(self, vw, ev, gp, psi_l, lai, cw, dt):
	    """Stored water flux, per unit ground area (um/s)""" 
	    return (vw - self.vwf(vw, ev, gp, psi_l, lai, cw, dt))*lai*10.**6/dt
	    #return self.gwf(self.psi_wf(vw,cw))*lai*(self.psi_wf(vw,cw) - ((ev*(1. - self.F_CAP))/(lai*gp) + psi_l))
	def qsf(self, soil, s, zr, psi_s, psi_b): ####fixed
	    """Soil water flux, per unit ground area (um/s)"""
	    return self.gsr(soil, s, zr, self.lam_s)*(psi_s-psi_b)
	def qgwf(self, soil, psi_gw, psi_b, gw_cs):############ Fixed
	    return self.ggwr(soil, soil.s, self.zr, self.lam_gw)*(psi_gw - psi_b)
	def qbx(self, gp, psi_x, psi_b, lai):
	    return (gp*lai/self.F_CAP)*(psi_b - psi_x)
	def gwf(self, psi_w):
	    """Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
	    return self.GWMAX*exp(-(-psi_w/2.)**2.)
	    #return GWMAX[species]*(vw/VWT[species])**4. 
	def gsrfp(self, soil, s, gp, lai, zr): #fix this delete!!!
	    """Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
	    #return (lai*self.gsr(soil, s, zr)*gp/self.F_CAP)/(self.gsr(soil, s, zr) +  lai*gp/self.F_CAP)
	def ggwrfp(self, soil, s, gp, lai, zr): #### Not sure about this - should LAI factor into groundwater paramete Delete this?
	    """Groundwater-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
	    #return (lai*self.ggwr(soil, s, zr)*gp/self.F_CAP)/(self.gsr(soil, s, zr) +  lai*gp/self.F_CAP)
	def vwf(self, vw, ev, gp, psi_l, lai, cw, dt): ##############lai on bottom of psi_ equation
	    """Stored water volume, per unit leaf area (m3/m2)"""
	    psi_w = self.psi_wf(vw, cw)
	    return min(vw - self.gwf(psi_w)*(psi_w - (ev*(1. - self.F_CAP))/(lai*gp) - psi_l)*dt/10.**6, self.VWT)
	def fBal(self, params, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, psi_w, gw_cs): ############Fix this
	    psi_l, tl = params
	    if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
		    return (phi*lai - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)-self.qbx(gp,self.psi_x(self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,gp,lai),self.psi_b(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,psi_w,gp,self.gwf(psi_w),lai,self.cw,soil),lai) - self.qwf(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),gp,psi_l,lai, self.cw,dt))\
			    #self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
    			#-((self.a(self.lai, soil.psi_s(s), soil.psi_gw(soil.gw_cs), self.gwf(self.psi_wf(self.vw, self.cw)), self.psi_wf(self.vw, self.cw), soil, soil.s, self.zr, self.lam_s, self.lam_gw) + \
    			#(self.d(soil, soil.s, self.zr) + self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)*self.psi_wf(self.vw, self.cw) - \
				#(self.d(soil, soil.s, self.zr)*(1 + self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)) + self.gwf(self.psi_wf(self.vw, self.cw))*self.lai)*psi_l)/ \
				#(1 + (self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)*((1/self.F_CAP) - 1)) + self.d(soil, soil.s, self.zr)*self.b(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)))))
	    else:
	    	return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) - self.qbx(gp,self.psi_x(self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,gp,lai), self.psi_b(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),psi_l,psi_w,gp,self.gwf(psi_w),lai,self.cw,soil),lai) - self.qwf(self.vw,self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared),gp,psi_l,lai, self.cw,dt))
		    	#self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) \
				#-((self.a(self.lai, soil.psi_s(s), soil.psi_gw(soil.gw_cs), self.gwf(self.psi_wf(self.vw, self.cw)), self.psi_wf(self.vw, self.cw), soil, soil.s, self.zr, self.lam_s, self.lam_gw) + \
    			#(self.d(soil, soil.s, self.zr) + self.e(self.gp, self.lai))*self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai)*self.psi_wf(self.vw, self.cw) - \
				#(((self.d(soil, soil.s, self.zr) + self.e(self.gp, self.lai))*(1 + (self.c(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai))) - self.e(self.gp, self.lai))*psi_l))/ \
				#((-(1-self.F_CAP)/self.F_CAP) + (self.d(soil, soil.s, self.zr) + self.e(self.gp, self.lai))*self.b(self.gp, self.gwf(self.psi_wf(self.vw, self.cw)), self.lai))))