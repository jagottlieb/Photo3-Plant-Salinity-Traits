from math import exp, pi, sqrt, log
from scipy.optimize import fsolve, least_squares
from sympy import *
import numpy as np
from dics import *
from functions import *

class Hydro(object):
	A_ROOT = 8.
	def __init__(self, species):
		self.GPMAX = species.GPMAX
		self.GA = species.GA
		self.gp = species.GPMAX
		self.GCUT = species.GCUT
		self.RAIW = species.RAIW
		self.zr = species.ZR
		self.lai = species.LAI
		self.psi_l_a = []
		self.gp_a = []
		self.gsv_a = []
		self.tl_a = []
		self.ev_a = []
	def rai(self, s):
		"""Root area index (-)"""
		return self.RAIW*s**-self.A_ROOT
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
	def gpf(self, psi_l):
		"""Plant conductance, per unit leaf area (um/(s-MPa))"""
		# if psi_l<-10:
		# 	return 0.
		# else:
		# 	return self.GPMAX*exp(-(-psi_l/2.)**2.)
		# return self.GPMAX
		# Add minimum threshold to prevent division by zero in downstream calculations
		gp_val = self.GPMAX*exp(-(-psi_l/2.)**2.)
		return max(gp_val, 1e-10)  # Ensure gp never goes below 1e-10
	def shf(self, tl, ta, lai):
		"""Sensible heat flux (W/m^2), per unit ground area"""
		return CP_A*RHO_A*self.GA*(tl-ta)/1000.*lai
	def gsr(self, soil, s, zr):
		"""Soil-Root Conductance, per unit ground area (um/(s-MPa))"""
		return (soil.leak(s)*sqrt(self.rai(s))*1000000.)/(float(pi)*g*RHO_W*zr)
	def gsw(self, photo, phi, ta, psi_l, qa, tl, ci, ared): 
		"""Stomatal conductance to water, per unit leaf area (mol/m2/sec)"""
		#return gsN(phi, Ta, psi_l, qa, Tl, ci, t)*(1.6*(1. + GMGSRATIO[pType[species]]))/(1.6 + GMGSRATIO[pType[species]]) + (gcut[species]*po/(1000.*R*Ta))
		return photo.gsc(phi, ta, psi_l, qa, tl, ci, ared)*1.6 + (self.GCUT*P_ATM/(1000.*R*ta))

class HydroNC(Hydro):
	def __init__(self, species, atm, soil, photo, vwi):
		Hydro.__init__(self, species)
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
# 		if self.qi(self.tl, self.psi_l) < atm.qa:
# 			self.psi_l = psi_i(atm.ta, atm.qa)
# 			self.tl = fsolve(self.fBal_psil_known, (290.), args= (self.psi_l, soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
	def update(self, atm, soil, photo, dt):
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
# 		if self.qi(self.tl, self.psi_l) < atm.qa:
# 			self.psi_l = psi_i(atm.ta, atm.qa)
# 			self.tl = fsolve(self.fBal_psil_known, (290.), args= (self.psi_l, soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.qs = self.ev
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
	def output(self):
		return {'psi_l': self.psi_l_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 'ev': self.ev_a}
	def gsrp(self, soil, s, gp, lai, zr):
		"""Soil-Root-Plant Conductance, per unit ground area (um/(s-MPa))"""
		return (lai*self.gsr(soil, s, zr)*gp)/(self.gsr(soil, s, zr) + lai*gp)
	def qsf(self, photo, phi, ta, psi_l, qa, tl, ci, lai):
		return self.evf(photo, phi, ta, psi_l, qa, tl, ci, lai, photo.ared)
	def fBal(self, p, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr):
		psi_l, tl =p

		if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
			return (phi*lai - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., \
				self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) - self.gsrp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l))
		else:
			return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., \
				self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) - self.gsrp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l))
	def fBal_psil_known(self, p, psi_l, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr):
		tl = p

		if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
			return phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.gsrp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l)/1000000.
		else:
			return phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.gsrp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l)/1000000.


class HydroCap(Hydro):
	F_CAP = 0.5


# # Below is the protocol for solving with Penman-Monteith (no solving energy balance!)
# 	def __init__(self, species, atm, soil, photo, vwi):
# 		Hydro.__init__(self, species)
# 		self.GWMAX = species.GWMAX
# 		self.VWT = species.VWT
# 		self.CAP = species.CAP
# 		self.vw = vwi*self.VWT
# 		self.tl = atm.ta
# 		self.psi_l = fsolve(self.fBalPen, (-1.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.tl))
# 		self.gp = self.gpf(self.psi_l)
# 		self.ev = self.evfPen(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
# 		self.tl = self.tlPen(atm.ta, atm.phi, self.ev)
# 		self.vw_a = []

# 	def update(self, atm, soil, photo, dt):
# 		self.psi_l = fsolve(self.fBalPen, (-1.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.tl))
# 		self.gp = self.gpf(self.psi_l)
# 		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, dt)
# 		self.ev = self.evfPen(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
# 		self.tl = self.tlPen(atm.ta, atm.phi, self.ev)
# 		self.qs = self.qsf(self.vw, self.ev, self.gp, self.psi_l, self.lai, dt)
# 		self.psi_l_a.append(self.psi_l)
# 		self.gp_a.append(self.gp)
# 		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
# 		self.tl_a.append(self.tl) 
# 		self.ev_a.append(self.ev)
# 		self.vw_a.append(self.vw)

	def __init__(self, species, atm, soil, photo, vwi):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.CAP = species.CAP
		self.vw = vwi*self.VWT
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
# 		if self.qi(self.tl, self.psi_l) < atm.qa:
# 			self.psi_l = psi_i(atm.ta, atm.qa)
# 			self.tl = fsolve(self.fBal_psil_known, (290.), args= (self.psi_l, soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.vw_a = []

	def update(self, atm, soil, photo, dt):
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
# 		if self.qi(self.tl, self.psi_l) < atm.qa:
# 			self.psi_l = psi_i(atm.ta, atm.qa)
# 			self.tl = fsolve(self.fBal_psil_known, (290.), args= (self.psi_l, soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
		self.gp = self.gpf(self.psi_l)
		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, dt)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.qs = self.qsf(self.vw, self.ev, self.gp, self.psi_l, self.lai, dt)
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
		self.vw_a.append(self.vw)

	def output(self):
		return {'psi_l': self.psi_l_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 'ev': self.ev_a, 'vw': self.vw_a}

	def psi_wf(self, vw): 
		"""Water potential of stored water (MPa)"""
		return (1./self.CAP)*vw/self.VWT - (1./self.CAP)
	def vwf(self, vw, ev, gp, psi_l, lai, dt):
		"""Stored water volume, per unit leaf area (m3/m2)"""
		if gp < 0.000001:
			return vw
		else:
			return min(vw - self.gwf(self.psi_wf(vw))*(self.psi_wf(vw) - (ev*(1. - self.F_CAP))/(lai*gp) - psi_l)*dt/10.**6, self.VWT)
	def psi_xf(self, ev, gp, psi_l):
		"""Water potential at connection node x (MPa)"""
		return ev*(1. - self.F_CAP)/(lai*gp) + psi_l
	def qwf(self, vw, ev, gp, psi_l, lai, dt):
		"""Stored water flux, per unit ground area"""
		return (vw - self.vwf(vw, ev, gp, psi_l, lai, dt))*lai*10.**6/dt
	def qsf(self, vw, ev, gp, psi_l, lai, dt):
		"""Soil water flux, per unit ground area"""
		return ev - self.qwf(vw, ev, gp, psi_l, lai, dt)
	def gwf(self, psi_w):
		"""Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
		# return self.GWMAX
		return self.GWMAX*exp(-(-psi_w/2.)**2.)
		#return GWMAX[species]*(vw/VWT[species])**4. 
	def gsrfp(self, soil, s, gp, lai, zr):
		"""Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
		return (lai*self.gsr(soil, s, zr)*gp/self.F_CAP)/(self.gsr(soil, s, zr) +  lai*gp/self.F_CAP)
	def fBal(self, params, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw):
		psi_l, tl = params
		psi_w = self.psi_wf(vw)

		if gp == 0.:
			return (psi_l - psi_i(ta, qa), \
				phi*lai - self.shf(tl, ta, lai))
		elif lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
			return (phi*lai - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000.,  \
				self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
				-(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))
		else:
			# energy balance, phi = shf + lambda rho evf
			# water balance, evf = (gsrfp(psi_s - psi_l) + gw(psi_w-psi_l))/(1 + etc....)
			return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., \
				self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
				-(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))
	def fBal_psil_known(self, p, psi_l, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw):
		tl = p
		psi_w = self.psi_wf(vw)
		if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
			return phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp)/1000000.
		else:
			return phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp)/1000000.
	def tlPen(self, ta, phi, ev):
		H = phi - LAMBDA_W*RHO_W*ev/1000000.
		return ta + H*1000/(CP_A*RHO_A*self.GA)

	def fBalPen(self, param, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw, tl):
		psi_l = param
		psi_w = self.psi_wf(vw)
		return self.evfPen(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) \
		- (self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
		(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp)

# Halophyte class copied from GWHalophyteGW, renamed for multi-compartment support
class Halophyte(Hydro):
	F_CAP = 0.5
	E = 0.90 # filtration efficiency, unitless
	Salt_Uptake = True
	TS = 293. # soil water temp (K)
	GWTS = 293 # groundwater temp (K)
	IV = 2. # van't hoff coefficient for NaCl
	#psi_w_i = -1.8
	TL = 293

	def __init__(self, species, atm, soil, photo, vwi, cw, s_arr, root_frac_arr, B, cs_arr, wr, wft, pi0, eta, mcap, dt):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.vw = vwi * self.VWT
		self.MW = cw * self.vw
		self.CAP = species.CAP
		self.cw = self.MW / self.vw
		self.delta_cw = 0
		self.delta_psi_w = 0
		self.w = self.vw / self.VWT
		self.psi_w_0 = self.cw * R * self.TS / (self.IV * self.TL * 10 ** (-6.))
		self.hr_cum = 0
		self.wr = wr
		self.wft = wft
		self.pi0 = pi0
		self.eta = eta
		self.mcap = mcap
		self.dt = dt

		# Arrays for time series outputs
		self.vw_a = []
		self.cw_a = []
		self.psi_b_a = []
		self.psi_x_a = []
		self.qw_a = []
		self.w_a = []
		self.psi_w_a = []
		self.psi_w_osm_a = []
		self.psi_w_turgor_a = []
		self.gs_a = []
		self.gwf_a = []
		self.gsr_a = [[] for _ in s_arr]
		self.qs_a = []
		self.hr_cum_a = []
		self.energy_balance = []
		self.water_balance = []
		self.flux_balance = []
		self.uptake_a = []
		self.MW_a = []
		self.delta_psi_w_a = []
		self.ev_a = []
		self.psi_l_a = []
		self.tl_a = []
		self.gp_a = []
		self.gsw_a = []
		self.qbx_a = []

		# Initial solve for psi_l and tl
		self.out = least_squares(
			self.fBal,
			(-1, 292.),
			args=(
				soil, photo, atm.phi, atm.ta, atm.qa, photo.cm,
				s_arr, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw),
				root_frac_arr, B, cs_arr, dt
			),
			bounds=([-10.0, 260.0], [0.0, 330.0]),
			method='trf', ftol=3e-16, xtol=3e-16, x_scale='jac', max_nfev=2000
		)
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)

	def update(self, atm, soil, photo, root_frac_arr, B, cs_arr, dt):
		self.cw = self.MW / self.vw
		self.out = least_squares(
			self.fBal,
			(-1, 292.),
			args=(
				soil, photo, atm.phi, atm.ta, atm.qa, photo.cm,
				soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw),
				root_frac_arr, B, cs_arr, dt
			),
			bounds=([-10.0, 260.0], [0.0, 330.0]),
			method='trf', ftol=3e-16, xtol=3e-16, x_scale='jac', max_nfev=2000
		)
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)
		psi_b_val = self.psi_b(soil, soil.s, self.zr, self.psi_l, soil.psi_s(soil.s,soil.cs), B, root_frac_arr, self.gp, self.lai, self.ev, self.vw, self.cw)
		self.qs = self.qsf(soil, soil.s, self.zr, soil.psi_s(soil.s,soil.cs), psi_b_val, B, root_frac_arr)
		self.gw = self.gwf(self.psi_wf(self.vw, self.cw))
		self.qw = self.qwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)
		self.qbxf = self.qbx(
			self.gp,
			self.psi_x(self.ev, self.psi_l, self.gp, self.lai),
			self.psi_b(soil, soil.s, self.zr, self.psi_l, soil.psi_s(soil.s,soil.cs), B, root_frac_arr, self.gp, self.lai, self.ev, self.vw, self.cw),
			self.lai)
		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)
		self.w = self.vw / self.VWT
		self.flux_balance.append(np.sum(self.qs)+self.qwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)-self.ev)
		self.uptake_val = self.uptake(self.qsf(soil, soil.s, self.zr, soil.psi_s(soil.s,soil.cs), psi_b_val, B, root_frac_arr), cs_arr)
		self.hr_cum = self.hr_cum + np.sum(self.hr(self.qs)) * 30 * 60 / 1000
		# Store time series outputs
		self.qs_a.append(self.qs)
		gsr_vals = self.gsr(soil, soil.s, self.zr, B, root_frac_arr)
		for i, gsr_val in enumerate(gsr_vals):
			self.gsr_a[i].append(gsr_val)
		self.psi_l_a.append(self.psi_l)
		self.tl_a.append(self.tl)
		self.ev_a.append(self.ev)
		self.qw_a.append(self.qw)
		self.qbx_a.append(self.qbxf)
		self.gp_a.append(self.gp)
		self.gwf_a.append(self.gw)
		self.gsw_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
		self.vw_a.append(self.vw)
		self.cw_a.append(self.cw)
		#self.psi_w_a.append(self.psi_wf( self.vw, self.cw))
		self.psi_w_a.append(
			self.psi_wf_bartlett(
				self.vw, VWT=self.VWT, pi_0=self.pi0, wft=self.wft, wr=self.wr, eta=self.eta, psi_0=0, mcap=self.mcap
			)
			)
		self.psi_w_turgor_a.append(
			self.psi_wf_turgor_bartlett(
				self.vw, VWT=self.VWT, pi_0=self.pi0, wft=self.wft, wr=self.wr, eta=self.eta,
				)
				)
		self.psi_w_osm_a.append(
			self.psi_wf_osm_bartlett(
				self.vw, VWT=self.VWT, pi_0=self.pi0, wft=self.wft, wr=self.wr,
				)
				)
		self.psi_x_a.append(self.psi_x(self.ev, self.psi_l, self.gp, self.lai))
		self.psi_b_a.append(self.psi_b(soil, soil.s, self.zr, self.psi_l, soil.psi_s(soil.s,soil.cs), B, root_frac_arr, self.gp, self.lai, self.ev, self.vw, self.cw))
		self.hr_cum_a.append(self.hr_cum)
		self.w_a.append(self.w)
		self.MW_a.append(self.MW)
		# self.flux_balance.append(np.sum(self.qs)+self.qwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)-self.ev)

	def output(self):
		return {
			'psi_l': self.psi_l_a, 
			'qs': self.qs_a,
			'gp': self.gp_a,
			'gsr': self.gsr_a,
			'gsw': self.gsw_a, 
			'tl': self.tl_a, 
			'ev': self.ev_a, 
			'ev_cum': np.cumsum(list(i*1.8 for i in self.ev_a)),
			'qw': self.qw_a,
			'vw': self.vw_a, 
			'cw': self.cw_a, 
			'psi_w': self.psi_w_a,
			'psi_w_turgor': self.psi_w_turgor_a,
			'psi_w_osm': self.psi_w_osm_a,
			'psi_x': self.psi_x_a,
			'psi_b': self.psi_b_a,
			'gwf': self.gwf_a,
			'w': self.w_a,
			'flux_balance': self.flux_balance,
			'hr_cum': self.hr_cum_a, 
			'Uptake': self.uptake_a, 
			'MW': self.MW_a, 
			'delta psi w': self.delta_psi_w_a,
			'qbx': self.qbx_a,
		}
	def psi_wf_bartlett(
			self,
			vw, #Volumetric water content in sapwood
			VWT, # Volumetric water content in sapwood at saturation
			pi_0=-1.5, # Osmotic potential at full turgor (MPa)
        	wft=1, # Relative water content at full turgor
        	wr=0.1, # Relative water content at point of apoplastic storage only
        	eta=5, # Bulk elastic modulus. Bartlett 2012 for a crop
        	psi_0=0, # Plant water storage potential - assumed to be 0 at full saturation. Bartlett 2012 for Mangroves
        	mcap=12 # Slope of PV curve between saturation and full turgor, from Bartlett et al., 2012
	):
		"""Calculate water potential using Bartlett et al. 2012 framework for a single value."""
		
		# Calculate key thresholds
		vr = VWT * wr  # Residual volumetric water content
		wtlp_tot = (1 - (vr/VWT)) * ((pi_0 + eta) / eta) + (vr/VWT)  # From Bartlett et al., 2012
		vtlp = wtlp_tot * VWT  # Volumetric water content at turgor loss point
		vft = wft * VWT  # Volumetric water content at full turgor
		# Clamp to physical range to avoid NaNs in downstream balance equations
		vw_safe = min(max(vw, vr + 1e-9), VWT)
		
		# Determine which condition applies and calculate psi_x
		if vft < vw_safe <= VWT:
			# First case: ψx = ψ0,x + mcap * (θx / θsat,x - 1)
			psi_wf = psi_0 + mcap * ((vw_safe / VWT) - 1)
		
		elif vtlp < vw_safe <= vft:
			# Second case: osmotic + turgor pressure
			psi_wf_osm = self.psi_wf_osm_bartlett(vw_safe, VWT, pi_0, wft, wr)
			psi_wf_turgor = self.psi_wf_turgor_bartlett(vw_safe, VWT, pi_0, wft, wr, eta)
			psi_wf = psi_wf_osm + psi_wf_turgor
		
		elif vr < vw_safe <= vtlp:
			# Third case: osmotic potential only (no turgor)
			psi_wf_osm = self.psi_wf_osm_bartlett(vw_safe, VWT, pi_0, wft, wr)
			psi_wf = psi_wf_osm
		
		else:
			# Out of range - fall back to clamped value
			psi_wf = self.psi_wf_osm_bartlett(vw_safe, VWT, pi_0, wft, wr)
		
		return psi_wf
	def psi_wf(self, vw, cw): 
		TL = 293.
		return self.psi_w_0 - (self.psi_w_0/self.w**(1/400)) - self.delta_cw*R*self.IV*TL*10.**(-6.) #+ 0.5*Plant_h*g*RHO_W*10**(-6)
	
	def psi_wf_osm_bartlett(self, vw, VWT, pi_0, wft, wr):
		"""Helper: Calculate osmotic potential using Bartlett et al. 2012 framework."""
		vr = VWT * wr  # Residual volumetric water content
		vw_safe = max(vw, vr + 1e-9)
		psi_osm = pi_0 * (VWT * wft - vr) / (vw_safe - vr)
		return psi_osm
	
	def psi_wf_turgor_bartlett(self, vw, VWT, pi_0, wft, wr, eta):
		"""Helper: Calculate turgor pressure using Bartlett et al. 2012 framework."""
		vr = VWT * wr  # Residual volumetric water content
		psi_turgor = abs(pi_0) - eta * (VWT * wft - vw) / (VWT * wft - vr)
		return max(psi_turgor, 0)  # Turgor cannot be negative
	
	def psi_wf_turgor(self, vw, cw=None, VWT=None, eta=27.7, aF=0.75, TL=293.):
		"""Calculate turgor pressure (MPa) based on salt concentration model."""
		if VWT is None:
			VWT = self.VWT
		psi_wf_osm_ft = -(self.MW/VWT)*R*self.IV*TL*10.**(-6.)
		return max(-psi_wf_osm_ft - eta*((1-(vw/VWT))/(1-aF)), 0)
	
	def psi_wf_osm(self, vw, cw=None, VWT=None, TL=293.):
		"""Calculate osmotic potential (MPa) based on salt concentration model."""
		if VWT is None:
			VWT = self.VWT
		psi_wf_osm_ft = -(self.MW/VWT)*R*self.IV*TL*10.**(-6.)
		return psi_wf_osm_ft/(vw/VWT)
	def a(self, soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr):
		"""Weighted sum of soil-root conductance times soil water potential for all compartments."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return np.sum(gsr_vals * np.array(psi_s_arr))
	def b(self, gp, gw, lai):
		return ((self.F_CAP/(gp*lai)) + ((1-self.F_CAP)/(gp*lai)) + ((self.F_CAP*(1-self.F_CAP)*gw/((gp**2)*lai))))
	def c(self, gp, gw, lai):
		return (self.F_CAP*gw/gp)
	def d(self, soil, s_arr, zr, psi_l, B, root_frac_arr):
		"""Sum of soil-root conductances for all compartments."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return np.sum(gsr_vals)
	def e(self, gp, lai):
		return (gp*lai/self.F_CAP)
	def psi_x(self, ev, psi_l, gp, lai): 
		# Add safeguard for very small gp values
		gp_safe = max(gp, 1e-10)
		return (ev*(1-self.F_CAP)/(lai*gp_safe) + psi_l)
	def psi_b(self, soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr, gp, lai, ev, vw, cw):
		"""Root/soil interface potential: algebraically enforces sum(qs_i) = qbx."""
		psi_x_val = self.psi_x(ev, psi_l, gp, lai)
		a_val = self.a(soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr)
		d_val = self.d(soil, s_arr, zr, psi_l, B, root_frac_arr)

		gp_term = gp * lai / self.F_CAP
		denom = d_val + gp_term
		if abs(denom) < 1e-12:
			denom = 1e-12 if denom >= 0 else -1e-12
		return (a_val + gp_term * psi_x_val) / denom
	def qwf(self, vw, ev, gp, psi_l, lai, cw, dt):
		return (vw - self.vwf(vw, ev, gp, psi_l, lai, cw, dt))*lai*10.**6/dt
	def qsf(self, soil, s_arr, zr, psi_s_arr, psi_b, B, root_frac_arr):
		"""Soil water flux for multiple compartments (array output)"""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return gsr_vals * (np.array(psi_s_arr) - psi_b)
	def hr(self, qsf_arr):
		return np.where(qsf_arr < 0, -qsf_arr, 0)
	def qbx(self, gp, psi_x, psi_b, lai):
		return (gp*lai/self.F_CAP)*(psi_b - psi_x)
	def gwf(self, psi_w):
		return self.GWMAX*exp(-(-psi_w/2.)**2.)
		#return self.GWMAX*(self.vw/self.VWT)**4
	def gsr(self, soil, s_arr, zr, B, root_frac_arr):
		"""Soil-Root Conductance for multiple compartments (array output, B is constant)"""
		rr = 0.2 * 10 ** -3
		kr = 10 ** -8
		gsr_list = []
		for s, root_frac in zip(s_arr, root_frac_arr):
			B_val = B * root_frac
			Ar = 2 * float(pi) * rr * B_val
			l = 0.53 / (float(pi) * B_val) ** 0.5
			ks = soil.leak(s) * 10 ** -6 / l
			gsr_val = (kr * ks / (kr + ks)) * 101.9 * 10 ** 6 * Ar * zr
			gsr_list.append(gsr_val)
		return np.array(gsr_list)
	def vwf(self, vw, ev, gp, psi_l, lai, cw, dt): 
		#psi_w = self.psi_wf(vw, cw)
		psi_w = self.psi_wf_bartlett(
			vw, VWT=self.VWT, pi_0=self.pi0, wft=self.wft, wr=self.wr, eta=self.eta, psi_0=0, mcap=self.mcap
		)
		# Add safeguard for very small gp values
		gp_safe = max(gp, 1e-10)
		return min(vw - self.gwf(psi_w)*(psi_w - (ev*(1. - self.F_CAP))/(lai*gp_safe) - psi_l)*dt/10.**6, self.VWT)
	def uptake(self, qsf_arr, cs_arr, E=None):
		"""Total salt uptake as the sum of qsf * cs across compartments, with conversion factors."""
		if E is None:
			E = self.E
		if self.Salt_Uptake:
			return np.sum(np.array(qsf_arr) * np.array(cs_arr)) * (1-E) * 30 * 60 * 10**(-6)
		else:
			return 0
	def fBal(
		self, params, soil, photo, phi, ta, qa, c1, s_arr, lai, gp, ared, zr, psi_w, root_frac_arr, B, cs_arr, dt
	):
		"""Energy and water balance equations for multi-compartment Halophyte."""
		psi_l, tl = params
		# Calculate transpiration
		evf_val = self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)
		# Calculate psi_s for each compartment
		psi_s_arr = soil.psi_s(soil.s, soil.cs)
		# Calculate psi_b (root/soil interface potential)
		psi_b_val = self.psi_b(
			soil, soil.s, self.zr, psi_l, soil.psi_s(soil.s,soil.cs), B, root_frac_arr, gp, lai, evf_val, self.vw, self.cw
		)
		# Calculate qsf (soil water flux for each compartment)
		qsf_arr = self.qsf(
			soil, soil.s, self.zr, soil.psi_s(soil.s,soil.cs), psi_b_val, B, root_frac_arr
		)
		# Calculate xylem-basal flux
		psi_x_val = self.psi_x(evf_val, psi_l, gp, lai)
		qbx_val = self.qbx(gp, psi_x_val, psi_b_val, lai)

		# Calculate stored water flux (single value)
		qwf_val = self.qwf(self.vw, evf_val, gp, psi_l, lai, self.cw, dt)
		# Energy balance equation
		energy_balance = (
			phi * lai
			- self.shf(tl, ta, lai)
			- LAMBDA_W * RHO_W * evf_val / 1_000_000.0
		)
		# Water balance equation
		water_balance = (
			evf_val
			- qbx_val
			- qwf_val
		)
		return (water_balance, energy_balance)
