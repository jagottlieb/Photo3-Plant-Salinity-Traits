from math import exp, pi, sqrt, log
from scipy.optimize import fsolve, least_squares
from sympy import *
import numpy as np
from dics import *
from functions import *

class Hydro(object):
	A_ROOT = 8.
	def __init__(self, species, gcut=None):
		self.GPMAX = species.GPMAX
		self.GA = species.GA
		self.gp = species.GPMAX
		self.GCUT = species.GCUT if gcut is None else gcut
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

# Halophyte class with multi-soil compartment capability and stem water storage
class HalophyteStemStorageMultiComp(Hydro):
	TS = 293. # soil water temp (K)
	GWTS = 293 # groundwater temp (K)
	IV = 2. # van't hoff coefficient for NaCl
	#psi_w_i = -1.8
	TL = 293 #

	def __init__(self, species, atm, soil, photo, vwi_stem, c_stem, s_arr, root_frac_arr, B, cs_arr, wr_stem, wft_stem, pi0_stem, eta_stem, mcap_stem, dt, salt_uptake=False, psi_wf_mode='bartlett', E=0.99, F_CAP=0.5, dynamic_E=False, c_stem_max=None):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.vw = vwi_stem * self.VWT
		self.CAP = species.CAP
		self.la = getattr(species, 'LA', None)
		self.delta_c_stem = 0
		self.delta_psi_w = 0
		self.w = self.vw / self.VWT
		self.hr_cum = 0
		self.wr_stem = wr_stem
		self.wft_stem = wft_stem
		self.pi0_stem = pi0_stem
		# c_stem is tracked as concentration in storage water, mol/m3 on a ground area basis.
		self.c_stem = c_stem
		self.MW = self.c_stem * self.vw
		self.psi_w = self.c_stem * R * self.TS *self.IV * 10 ** (-6.)
		self.eta_stem = eta_stem
		self.mcap_stem = mcap_stem
		self.dt = dt
		self.Salt_Uptake = salt_uptake
		# Reversible selector for plant storage water potential formulation.
		# Supported: 'bartlett' (default), 'legacy'.
		self.psi_wf_mode = psi_wf_mode
		self.E = E
		self.F_CAP = F_CAP
		self.dynamic_E = dynamic_E
		# c_stem_max uses the same concentration basis as c_stem: mol/m3 on a ground area basis.
		self.c_stem_max = c_stem_max

		# Arrays for time series outputs
		self.vw_a = []
		self.c_stem_a = []
		self.psi_b_a = []
		self.psi_x_a = []
		self.qw_stem_a = []
		self.w_stem_a = []
		self.psi_w_stem_a = []
		self.psi_w_osm_a = []
		self.psi_w_turgor_a = []
		self.gs_a = []
		self.gw_stem_a = []
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
		self.qw_a = self.qw_stem_a
		self.w_a = self.w_stem_a
		self.psi_w_a = self.psi_w_stem_a
		self.gwf_a = self.gw_stem_a

		# Initial solve for psi_l and tl
		self.out = least_squares(
			self.fBal,
			(-1, 292.),
			# args=(
			# 	soil, photo, atm.phi, atm.ta, atm.qa, photo.cm,
			# 	s_arr, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.c_stem),
			# 	root_frac_arr, B, cs_arr, dt
			# ),
			args=(
				soil, photo, atm.phi, atm.ta, atm.qa, photo.cm,
				s_arr, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.c_stem),
				root_frac_arr, B, cs_arr, dt
			),
			bounds=([-10.0, 260.0], [0.0, 330.0]),
			method='trf', ftol=3e-16, xtol=3e-16, x_scale='jac', max_nfev=2000
		)
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)

	def update(self, atm, soil, photo, root_frac_arr, B, cs_arr, dt):
		self.out = least_squares(
			self.fBal,
			(-1, 292.),
			args=(
				soil, photo, atm.phi, atm.ta, atm.qa, photo.cm,
				soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.c_stem),
				root_frac_arr, B, cs_arr, dt
			),
			bounds=([-10.0, 260.0], [0.0, 330.0]),
			method='trf', ftol=3e-16, xtol=3e-16, x_scale='jac', max_nfev=2000
		)
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)
		psi_b_val = self.psi_b(soil, soil.s, self.zr, self.psi_l, soil.psi_s(soil.s,soil.cs), B, root_frac_arr, self.gp, self.lai, self.ev, self.vw, self.c_stem)
		self.qs = self.qsf(soil, soil.s, self.zr, soil.psi_s(soil.s,soil.cs), psi_b_val, B, root_frac_arr)
		self.gw_stem = self.gwf_stem(self.psi_wf_stem(self.vw, self.c_stem))
		self.qw_stem = self.qwf_stem(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.c_stem, dt)
		# Backward-compatible scalar aliases.
		self.gw = self.gw_stem
		self.qw = self.qw_stem
		self.qbxf = self.qbx(
			self.gp,
			self.psi_x(self.ev, self.psi_l, self.gp, self.lai),
			self.psi_b(soil, soil.s, self.zr, self.psi_l, soil.psi_s(soil.s,soil.cs), B, root_frac_arr, self.gp, self.lai, self.ev, self.vw, self.c_stem),
			self.lai)
		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.c_stem, dt)
		self.w = self.vw / self.VWT
		self.flux_balance.append(np.sum(self.qs)+self.qwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.c_stem, dt)-self.ev)
		self.uptake_val = self.uptake(self.qsf(soil, soil.s, self.zr, soil.psi_s(soil.s,soil.cs), psi_b_val, B, root_frac_arr), cs_arr)
		self.MW = self.MW + self.uptake_val * dt
		# Updated c_stem concentration, mol/m3 on a ground area basis.
		self.c_stem = self.MW / self.vw
		self.hr_cum = self.hr_cum + np.sum(self.hr(self.qs)) * 30 * 60 / 1000
		#=========================================
		# Store time series outputs
		#=========================================
		self.qs_a.append(self.qs)
		gsr_vals = self.gsr(soil, soil.s, self.zr, B, root_frac_arr)
		for i, gsr_val in enumerate(gsr_vals):
			self.gsr_a[i].append(gsr_val)
		self.psi_l_a.append(self.psi_l)
		self.tl_a.append(self.tl)
		self.ev_a.append(self.ev)
		self.qw_stem_a.append(self.qw_stem)
		self.qbx_a.append(self.qbxf)
		self.gp_a.append(self.gp)
		self.gw_stem_a.append(self.gw_stem)
		self.gsw_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
		self.vw_a.append(self.vw)
		self.c_stem_a.append(self.c_stem)
		#self.psi_w_stem_a.append(self.psi_wf_stem(self.vw, self.c_stem))
		self.psi_w_stem_a.append(self.psi_wf_stem(self.vw, self.c_stem))
		self.psi_w_turgor_a.append(
			self.psi_wf_turgor_bartlett(
				self.vw, VWT=self.VWT, pi_0=self.pi0_stem, wft=self.wft_stem, wr=self.wr_stem, eta=self.eta_stem,
				)
				)
		self.psi_w_osm_a.append(
			self.psi_wf_osm_bartlett(
				self.vw, VWT=self.VWT, pi_0=self.pi0_stem, wft=self.wft_stem, wr=self.wr_stem,
				)
				)
		self.psi_x_a.append(self.psi_x(self.ev, self.psi_l, self.gp, self.lai))
		self.psi_b_a.append(self.psi_b(soil, soil.s, self.zr, self.psi_l, soil.psi_s(soil.s,soil.cs), B, root_frac_arr, self.gp, self.lai, self.ev, self.vw, self.c_stem))
		self.hr_cum_a.append(self.hr_cum)
		self.w_stem_a.append(self.w)
		self.MW_a.append(self.MW)
		# self.flux_balance.append(np.sum(self.qs)+self.qwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.c_stem, dt)-self.ev)

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
			'qw_stem': self.qw_stem_a,
			'vw': self.vw_a, 
			'c_stem': self.c_stem_a,
			'psi_w_stem': self.psi_w_stem_a,
			'psi_w_turgor': self.psi_w_turgor_a,
			'psi_w_osm': self.psi_w_osm_a,
			'psi_x': self.psi_x_a,
			'psi_b': self.psi_b_a,
			'gw_stem': self.gw_stem_a,
			'w_stem': self.w_stem_a,
			'flux_balance': self.flux_balance,
			'hr_cum': self.hr_cum_a, 
			'Uptake': self.uptake_a, 
			'MW': self.MW_a, 
			'delta psi w': self.delta_psi_w_a,
			'qbx': self.qbx_a,
			# Backward-compatible aliases
			'qw': self.qw_stem_a,
			'w': self.w_stem_a,
			'psi_w': self.psi_w_stem_a,
			'gwf': self.gw_stem_a,
		}

	def set_psi_wf_mode(self, mode):
		"""Set storage water potential mode at runtime: 'bartlett' or 'legacy'."""
		mode_norm = str(mode).strip().lower()
		if mode_norm not in ('bartlett', 'legacy'):
			raise ValueError("psi_wf_mode must be 'bartlett' or 'legacy'")
		self.psi_wf_mode = mode_norm

	def psi_wf_stem(self, vw, c_stem=None):
		"""Stem-suffixed wrapper for storage water potential."""
		return self.psi_wf(vw, c_stem=c_stem)

	def psi_wf(self, vw, c_stem=None):
		"""Output storage water potential formulation based on configured mode."""
		mode = str(getattr(self, 'psi_wf_mode', 'bartlett')).strip().lower()
		if mode == 'legacy':
			# Legacy behavior based on osmotic + turgor components.
			return self.psi_wf_osm(vw, c_stem=c_stem, VWT=self.VWT, TL=self.TS) + self.psi_wf_turgor(vw, c_stem=c_stem, VWT=self.VWT, TL=self.TS)
		return self.psi_wf_bartlett(
			vw,
			VWT=self.VWT,
			pi_0=self.pi0_stem,
			wft=self.wft_stem,
			wr=self.wr_stem,
			eta=self.eta_stem,
			psi_0=0,
			mcap=self.mcap_stem,
		)

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
		vr = VWT * wr  # Residual (apoplastic)volumetric water content
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
	# def psi_wf(self, vw, c_stem): 
	# 	TL = 293.
	# 	return self.psi_w_0 - (self.psi_w_0/self.w**(1/400)) - self.delta_c_stem*R*self.IV*TL*10.**(-6.) #+ 0.5*Plant_h*g*RHO_W*10**(-6)
	
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
	
	def psi_wf_turgor(self, vw, c_stem=None, VWT=None, eta=27.7, aF=0.75, TL=293.):
		"""Calculate turgor pressure (MPa) based on salt concentration model."""
		if VWT is None:
			VWT = self.VWT
		psi_wf_osm_ft = -(self.MW/VWT)*R*self.IV*TL*10.**(-6.)
		return max(-psi_wf_osm_ft - eta*((1-(vw/VWT))/(1-aF)), 0)
	
	def psi_wf_osm(self, vw, c_stem=None, VWT=None, TL=293.):
		"""Calculate osmotic potential (MPa) based on salt concentration model."""
		if VWT is None:
			VWT = self.VWT
		psi_wf_osm_ft = -(self.MW/VWT)*R*self.IV*TL*10.**(-6.)
		return psi_wf_osm_ft/(vw/VWT)
	def a(self, soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr):
		"""Aggregate soil-driven term for basal node solve (um/s)."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return np.sum(gsr_vals * np.array(psi_s_arr))
	def b(self, gp, gw, lai):
		return ((self.F_CAP/(gp*lai)) + ((1-self.F_CAP)/(gp*lai)) + ((self.F_CAP*(1-self.F_CAP)*gw/((gp**2)*lai))))
	def c(self, gp, gw, lai):
		return (self.F_CAP*gw/gp)
	def d(self, soil, s_arr, zr, psi_l, B, root_frac_arr):
		"""Total soil-root conductance across compartments (um/(s-MPa))."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return np.sum(gsr_vals)
	def e(self, gp, lai):
		return (gp*lai/self.F_CAP)
	def psi_x(self, ev, psi_l, gp, lai): 
		"""Xylem node water potential from transpiration partitioning (MPa)."""
		# Add safeguard for very small gp values
		gp_safe = max(gp, 1e-10)
		return (ev*(1-self.F_CAP)/(lai*gp_safe) + psi_l)
	def psi_b(self, soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr, gp, lai, ev, vw, c_stem):
		"""Root-base interface potential enforcing flux closure (MPa)."""
		psi_x_val = self.psi_x(ev, psi_l, gp, lai)
		a_val = self.a(soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr)
		d_val = self.d(soil, s_arr, zr, psi_l, B, root_frac_arr)

		gp_term = gp * lai / self.F_CAP
		denom = d_val + gp_term
		if abs(denom) < 1e-12:
			denom = 1e-12 if denom >= 0 else -1e-12
		return (a_val + gp_term * psi_x_val) / denom
	def qwf(self, vw, ev, gp, psi_l, lai, c_stem, dt):
		"""Stem storage flux on a ground area basis (um/s)."""
		return (vw - self.vwf(vw, ev, gp, psi_l, lai, c_stem, dt))*lai*10.**6/dt

	def qwf_stem(self, vw, ev, gp, psi_l, lai, c_stem, dt):
		"""Stem storage flux on a ground area basis (um/s)."""
		return self.qwf(vw, ev, gp, psi_l, lai, c_stem, dt)
	def qsf(self, soil, s_arr, zr, psi_s_arr, psi_b, B, root_frac_arr):
		"""Soil water flux by compartment from conductance and potential gradient (um/s)."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return gsr_vals * (np.array(psi_s_arr) - psi_b)
	def hr(self, qsf_arr):
		"""Hydraulic redistribution magnitude from negative soil fluxes (um/s)."""
		return np.where(qsf_arr < 0, -qsf_arr, 0)
	def qbx(self, gp, psi_x, psi_b, lai):
		"""Flux from root base to stem/xylem connection node on ground area basis (um/s)."""
		return (gp*lai/self.F_CAP)*(psi_b - psi_x)
	def gwf(self, psi_w):
		"""Stem storage-to-xylem conductance (um/(MPa s))."""
		#return self.GWMAX*exp(-(-psi_w/2.)**2.)
		return self.GWMAX*(self.vw/self.VWT)**4

	def gwf_stem(self, psi_w):
		"""Stem storage-to-xylem conductance (um/(MPa s))."""
		return self.gwf(psi_w)
	def gsr(self, soil, s_arr, zr, B, root_frac_arr):
		"""Soil-root conductance by compartment (um/(s-MPa))."""
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
	def vwf(self, vw, ev, gp, psi_l, lai, c_stem, dt): 
		"""Updated stem storage water depth/volume-equivalent state variable (m)."""
		psi_w = self.psi_wf(vw, c_stem)
		# Add safeguard for very small gp values
		gp_safe = max(gp, 1e-10)
		return min(vw - self.gwf(psi_w)*(psi_w - (ev*(1. - self.F_CAP))/(lai*gp_safe) - psi_l)*dt/10.**6, self.VWT)
	def uptake(self, qsf_arr, cs_arr, E=None):
		"""Salt uptake (mol/(m2 ground area s)); uses qsf (um/s), cs (mol/m3), and filtration efficiency E (-)."""
		if E is None:
			c_max = self.c_stem_max if self.c_stem_max is not None else max(np.max(cs_arr), 1e-6)
			if self.dynamic_E:
				c_thresh = self.E * c_max
				if self.c_stem < c_thresh:
					E = self.E
				else:
					E = min(self.E + (1 - self.E) * (self.c_stem - c_thresh) / (c_max - c_thresh), 0.9999)
			else:
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
			soil, soil.s, self.zr, psi_l, soil.psi_s(soil.s,soil.cs), B, root_frac_arr, gp, lai, evf_val, self.vw, self.c_stem
		)
		# Calculate qsf (soil water flux for each compartment)
		qsf_arr = self.qsf(
			soil, soil.s, self.zr, soil.psi_s(soil.s,soil.cs), psi_b_val, B, root_frac_arr
		)
		# Calculate xylem-basal flux
		psi_x_val = self.psi_x(evf_val, psi_l, gp, lai)
		qbx_val = self.qbx(gp, psi_x_val, psi_b_val, lai)

		# Calculate stored water flux (single value)
		qwf_val = self.qwf(self.vw, evf_val, gp, psi_l, lai, self.c_stem, dt)
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


class HalophyteStemLeafStorageMultiComp(Hydro):
	"""Multi-compartment halophyte hydraulics with both stem and leaf storage."""
	TS = 293.
	IV = 2.

	def __init__(
		self,
		species,
		atm,
		soil,
		photo,
		vwi_stem,
		c_stem,
		s_arr,
		root_frac_arr,
		B,
		cs_arr,
		wr_stem, 
		wft_stem,
		pi0_stem,
		eta_stem,
		mcap_stem,
		dt,
		salt_uptake=False,
		psi_wf_mode='bartlett',
		vwi_leaf=1.0,
		c_leaf=None,
		wr_leaf=0.46, # Reported by Rieger and Daniell for Pecan leaves
		wft_leaf=1.0,
		pi0_leaf=-1.5,
		eta_leaf=5,
		mcap_leaf=12,
		leaf_uptake_frac=0.5,
		gcut=None,
		E=0.99,
		F_CAP=0.5,
		dynamic_E=False,
		c_stem_max=None,
	):
		Hydro.__init__(self, species, gcut=gcut)
		self.E = E
		self.F_CAP = F_CAP
		self.dynamic_E = dynamic_E
		# c_stem_max uses the same concentration basis as c_stem: mol/m3 on a ground area basis.
		self.c_stem_max = c_stem_max
		# Stem storage parameters
		self.GWMAX_STEM = species.GWMAX
		self.VWTSTEM = species.VWT
		self.vw_stem = vwi_stem * self.VWTSTEM
		self.w_stem = self.vw_stem / self.VWTSTEM
		self.wr_stem = wr_stem
		self.wft_stem = wft_stem
		self.pi0_stem = pi0_stem
		self.pi0_stem_base = pi0_stem
		self.eta_stem = eta_stem
		self.mcap_stem = mcap_stem
		# c_stem is tracked as concentration in storage water, mol/m3 on a ground area basis.
		self.c_stem = c_stem

		# Leaf storage parameters
		self.la = getattr(species, 'LA', None)
		self.GWLEAF = getattr(species, 'GWLEAF', getattr(species, 'GWMAXLEAF', 0.0))
		self.VWTLEAF = getattr(species, 'VWTLEAF', 1e-6)
		self.vw_leaf = vwi_leaf * self.VWTLEAF
		self.w_leaf = self.vw_leaf / self.VWTLEAF
		self.wr_leaf = wr_leaf
		self.wft_leaf = wft_leaf
		self.pi0_leaf = pi0_leaf
		self.pi0_leaf_base = pi0_leaf
		self.eta_leaf = eta_leaf
		self.mcap_leaf = mcap_leaf
		if c_leaf is None:
			self.c_leaf = self.c_stem
		else:
			self.c_leaf = c_leaf
		# c_leaf is tracked as concentration in storage water, mol/m3 on a ground area basis.

		# Salt masses
		self.MW_stem = self.c_stem * self.vw_stem * self.la
		self.MW_leaf = self.c_leaf * self.vw_leaf * self.la

		self.dt = dt
		self.Salt_Uptake = salt_uptake
		self.psi_wf_mode = psi_wf_mode
		# Fraction of salt uptake allocated to the leaf vs stem on a volume weighted-basis
		self.leaf_uptake_frac = float(np.clip(leaf_uptake_frac*self.VWTLEAF/(self.VWTLEAF + self.VWTSTEM), 0.0, 1.0)) 

		# Cumulative uptake-tracked salt moles (mol m^-2 ground)
		self.MW_uptake_stem = 0.0
		self.MW_uptake_leaf = 0.0
		self.MW_uptake_total = 0.0
		self.psi_uptake_stem = 0.0
		self.psi_uptake_leaf = 0.0

		# Time series outputs
		self.qs_a = []
		self.gsr_a = [[] for _ in s_arr]
		self.psi_l_a = []
		self.tl_a = []
		self.ev_a = []
		self.gp_a = []
		self.gsw_a = []
		self.qbx_a = []
		self.qxl_a = []
		self.psi_x_a = []
		self.psi_b_a = []

		self.qw_stem_a = []
		self.gw_stem_a = []
		self.psi_w_stem_a = []
		self.psi_w_osm_stem_a = []
		self.psi_w_turgor_stem_a = []
		self.vw_stem_a = []
		self.c_stem_a = []
		self.pi0_stem_a = []
		self.psi_uptake_stem_a = []
		self.MW_uptake_stem_a = []

		self.qw_leaf_a = []
		self.gw_leaf_a = []
		self.psi_w_leaf_a = []
		self.psi_w_osm_leaf_a = []
		self.psi_w_turgor_leaf_a = []
		self.vw_leaf_a = []
		self.c_leaf_a = []
		self.pi0_leaf_a = []
		self.psi_uptake_leaf_a = []
		self.MW_uptake_leaf_a = []
		self.MW_uptake_total_a = []

		self.flux_balance = []

		# Initial solve for psi_l and tl
		self.out = least_squares(
			self.fBal,
			(-1., 292.),
			args=(
				soil, photo, atm.phi, atm.ta, atm.qa, photo.cm,
				s_arr, self.lai, self.gp, 1., self.zr,
				root_frac_arr, B, cs_arr, dt,
			),
			bounds=([-10.0, 260.0], [0.0, 330.0]),
			method='trf', ftol=3e-16, xtol=3e-16, x_scale='jac', max_nfev=2000
		)
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)

	def _psi_wf_bartlett(self, vw, VWT, pi_0, wft, wr, eta, mcap, psi_0=0):
		vr = VWT * wr
		wtlp_tot = (1 - (vr / VWT)) * ((pi_0 + eta) / eta) + (vr / VWT)
		vtlp = wtlp_tot * VWT
		vft = wft * VWT
		vw_safe = min(max(vw, vr + 1e-9), VWT)
		if vft < vw_safe <= VWT:
			return psi_0 + mcap * ((vw_safe / VWT) - 1)
		if vtlp < vw_safe <= vft:
			psi_osm = pi_0 * (VWT * wft - vr) / (vw_safe - vr)
			psi_turgor = max(abs(pi_0) - eta * (VWT * wft - vw_safe) / (VWT * wft - vr), 0)
			return psi_osm + psi_turgor
		return pi_0 * (VWT * wft - vr) / (vw_safe - vr)

	def _psi_wf_components_bartlett(self, vw, VWT, pi_0, wft, wr, eta, mcap, psi_0=0):
		"""Return osmotic and turgor components (MPa) for Bartlett PV formulation."""
		vr = VWT * wr
		wtlp_tot = (1 - (vr / VWT)) * ((pi_0 + eta) / eta) + (vr / VWT)
		vtlp = wtlp_tot * VWT
		vft = wft * VWT
		vw_safe = min(max(vw, vr + 1e-9), VWT)

		psi_osm = pi_0 * (VWT * wft - vr) / (vw_safe - vr)
		if vft < vw_safe <= VWT:
			psi_total = psi_0 + mcap * ((vw_safe / VWT) - 1)
			psi_turgor = max(psi_total - psi_osm, 0)
		elif vtlp < vw_safe <= vft:
			psi_turgor = max(abs(pi_0) - eta * (VWT * wft - vw_safe) / (VWT * wft - vr), 0)
		else:
			psi_turgor = 0.0
		return psi_osm, psi_turgor

	def _psi_wf_legacy(self, vw, MW, VWT, TL=293., eta=27.7, aF=0.75):
		psi_osm_ft = -(MW / VWT) * R * self.IV * TL * 10.**(-6.)
		psi_osm = psi_osm_ft / (vw / VWT)
		psi_turgor = max(-psi_osm_ft - eta * ((1 - (vw / VWT)) / (1 - aF)), 0)
		return psi_osm + psi_turgor

	def _psi_wf_components_legacy(self, vw, MW, VWT, TL=293., eta=27.7, aF=0.75):
		"""Return osmotic and turgor components (MPa) for legacy storage formulation."""
		psi_osm_ft = -(MW / VWT) * R * self.IV * TL * 10.**(-6.)
		psi_osm = psi_osm_ft / (vw / VWT)
		psi_turgor = max(-psi_osm_ft - eta * ((1 - (vw / VWT)) / (1 - aF)), 0)
		return psi_osm, psi_turgor

	def psi_components_stem(self, vw):
		mode = str(getattr(self, 'psi_wf_mode', 'bartlett')).strip().lower()
		if mode == 'legacy':
			return self._psi_wf_components_legacy(vw, self.MW_stem, self.VWTSTEM, TL=self.TS)
		return self._psi_wf_components_bartlett(
			vw, self.VWTSTEM, self.pi0_stem, self.wft_stem, self.wr_stem, self.eta_stem, self.mcap_stem
		)

	def psi_components_leaf(self, vw):
		mode = str(getattr(self, 'psi_wf_mode', 'bartlett')).strip().lower()
		if mode == 'legacy':
			return self._psi_wf_components_legacy(vw, self.MW_leaf, self.VWTLEAF, TL=self.TS)
		return self._psi_wf_components_bartlett(
			vw, self.VWTLEAF, self.pi0_leaf, self.wft_leaf, self.wr_leaf, self.eta_leaf, self.mcap_leaf
		)

	def psi_wf_stem(self, vw):
		mode = str(getattr(self, 'psi_wf_mode', 'bartlett')).strip().lower()
		if mode == 'legacy':
			return self._psi_wf_legacy(vw, self.MW_stem, self.VWTSTEM, TL=self.TS)
		return self._psi_wf_bartlett(vw, self.VWTSTEM, self.pi0_stem, self.wft_stem, self.wr_stem, self.eta_stem, self.mcap_stem)

	def psi_wf_leaf(self, vw):
		mode = str(getattr(self, 'psi_wf_mode', 'bartlett')).strip().lower()
		if mode == 'legacy':
			return self._psi_wf_legacy(vw, self.MW_leaf, self.VWTLEAF, TL=self.TS)
		return self._psi_wf_bartlett(vw, self.VWTLEAF, self.pi0_leaf, self.wft_leaf, self.wr_leaf, self.eta_leaf, self.mcap_leaf)

	def gwf_stem(self, psi_w_stem):
		"""Stem storage-to-xylem conductance (um/(MPa s))."""
		return self.GWMAX_STEM * (self.vw_stem / self.VWTSTEM)**4

	def gwf_leaf(self, psi_w_leaf):
		"""Leaf storage-to-xylem conductance (um/(MPa s))."""
		if self.VWTLEAF <= 0:
			return 0.0
		return self.GWLEAF * (self.vw_leaf / self.VWTLEAF)**4

	def vwf_stem(self, vw_stem, ev, gp, psi_l, lai, dt, psi_x=None):
		"""Updated stem storage state variable from storage flux closure (m)."""
		psi_w_stem = self.psi_wf_stem(vw_stem)
		if psi_x is None:
			psi_x = self.psi_x(ev, psi_l, gp, lai)
		return min(vw_stem - self.gwf_stem(psi_w_stem) * (psi_w_stem - psi_x) * dt / 10.**6, self.VWTSTEM)

	def qwf_stem(self, vw_stem, ev, gp, psi_l, lai, dt, psi_x=None):
		"""Stem storage flux on a ground area basis (um/s)."""
		return (vw_stem - self.vwf_stem(vw_stem, ev, gp, psi_l, lai, dt, psi_x=psi_x)) * lai * 10.**6 / dt
	# m3 water/m2 leaf area 
	def vwf_leaf(self, vw_leaf, psi_l, dt):
		"""Updated leaf storage state variable from storage flux closure (m)."""
		psi_w_leaf = self.psi_wf_leaf(vw_leaf)
		return min(vw_leaf - self.gwf_leaf(psi_w_leaf) * (psi_w_leaf - psi_l) * dt / 10.**6, self.VWTLEAF)

	def qwf_leaf(self, vw_leaf, psi_l, lai, dt):
		"""Leaf storage flux on a ground area basis (um/s)."""
		return (vw_leaf - self.vwf_leaf(vw_leaf, psi_l, dt)) * lai * 10.**6 / dt

	def psi_x(self, ev, psi_l, gp, lai, gw_leaf=None, psi_w_leaf=None):
		"""Xylem node water potential including leaf storage coupling (MPa)."""
		# Requested updated form including leaf storage coupling.
		gp_safe = max(gp, 1e-10)
		if gw_leaf is None:
			gw_leaf = self.gwf_leaf(self.psi_wf_leaf(self.vw_leaf))
		if psi_w_leaf is None:
			psi_w_leaf = self.psi_wf_leaf(self.vw_leaf)
		return ev * (1 - self.F_CAP) / (lai * gp_safe) + psi_l - gw_leaf * (psi_w_leaf - psi_l)

	def a(self, soil, s_arr, zr, psi_s_arr, B, root_frac_arr):
		"""Aggregate soil-driven term for basal node solve (um/s)."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return np.sum(gsr_vals * np.array(psi_s_arr))

	def d(self, soil, s_arr, zr, B, root_frac_arr):
		"""Total soil-root conductance across compartments (um/(s-MPa))."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return np.sum(gsr_vals)

	def psi_b(self, soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr, gp, lai, ev, gw_leaf, psi_w_leaf):
		"""Root-base interface potential enforcing flux closure (MPa)."""
		psi_x_val = self.psi_x(ev, psi_l, gp, lai, gw_leaf=gw_leaf, psi_w_leaf=psi_w_leaf)
		a_val = self.a(soil, s_arr, zr, psi_s_arr, B, root_frac_arr)
		d_val = self.d(soil, s_arr, zr, B, root_frac_arr)
		gp_term = gp * lai / self.F_CAP
		denom = d_val + gp_term
		if abs(denom) < 1e-12:
			denom = 1e-12 if denom >= 0 else -1e-12
		return (a_val + gp_term * psi_x_val) / denom

	# Soil water flux on a per-gound area basis 
	def qsf(self, soil, s_arr, zr, psi_s_arr, psi_b, B, root_frac_arr):
		"""Soil water flux by compartment from conductance and potential gradient (um/s)."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return gsr_vals * (np.array(psi_s_arr) - psi_b)

	def qbx(self, gp, psi_x, psi_b, lai):
		"""Flux from root base to stem/xylem connection node on ground area basis (um/s)."""
		return (gp * lai / self.F_CAP) * (psi_b - psi_x)

	def qxl(self, gp, psi_x, psi_l, lai):
		"""Xylem-to-leaf flux (um/s) for the stem+leaf storage configuration."""
		return gp * lai / (1 - self.F_CAP) * (psi_x - psi_l)

	def gsr(self, soil, s_arr, zr, B, root_frac_arr):
		"""Soil-root conductance by compartment (um/(s-MPa))."""
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

	def hr(self, qsf_arr):
		"""Hydraulic redistribution magnitude from negative soil fluxes (um/s)."""
		return np.where(qsf_arr < 0, -qsf_arr, 0)

	# Salt uptake (mol m^-2 ground area s^-1) calculated as sum of qsf * cs across compartments, with conversion factors and optional dynamic E adjustment.
	def uptake(self, qsf_arr, cs_arr, E=None):
		"""Salt uptake (mol/(m2 ground area s)); uses qsf (um/s), cs (mol/m3), and filtration efficiency E (-)."""
		if E is None:
			# c_stem and c_stem_max are treated as mol/m3 on a ground area basis in this E adjustment.
			c_max = self.c_stem_max if self.c_stem_max is not None else max(np.max(cs_arr), 1e-6)
			if self.dynamic_E:
				c_thresh = self.E * c_max
				if self.c_stem < c_thresh:
					E = self.E
				else:
					# E = min(self.E + (1 - self.E) * (self.c_stem - c_thresh) / (c_max - c_thresh), 0.9999)
					E = min(self.E + (1 - self.E) * (self.c_stem/c_max), 0.9999)
			else:
				E = self.E
		if self.Salt_Uptake:
			# qsf (um s^-1) * cs (mol m^-3) -> mol m^-2 s^-1 after 1e-6 factor
			# and positive root inflow only is treated as plant uptake.
			return np.sum(np.maximum(np.array(qsf_arr), 0.0) * np.array(cs_arr)) * (1 - E) * 10**(-6)
		return 0
	#Units: VWTSTEM: m^3 water/m^2 leaf, la: m^2 leaf area, lai: m^2 leaf area/m^2 ground
	def _storage_volume_stem(self):
		return max(self.VWTSTEM * self.lai, 1e-12)

	def _storage_volume_leaf(self):
		return max(self.VWTLEAF * self.lai, 1e-12)

	def _update_pi0_from_uptake(self):
		cs_uptake_stem = self.MW_uptake_stem / self._storage_volume_stem()
		cs_uptake_leaf = self.MW_uptake_leaf / self._storage_volume_leaf()

		# Van't Hoff osmotic contribution from cumulative uptake (negative MPa).
		self.psi_uptake_stem = -cs_uptake_stem * R * self.IV * self.TS * 10.**(-6.)
		self.psi_uptake_leaf = -cs_uptake_leaf * R * self.IV * self.TS * 10.**(-6.)

		# Effective pi0 becomes more negative as uptake accumulates.
		self.pi0_stem = self.pi0_stem_base + self.psi_uptake_stem
		self.pi0_leaf = self.pi0_leaf_base + self.psi_uptake_leaf

	def fBal(self, params, soil, photo, phi, ta, qa, c1, s_arr, lai, gp, ared, zr, root_frac_arr, B, cs_arr, dt):
		psi_l, tl = params
		evf_val = self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)

		psi_s_arr = soil.psi_s(soil.s, soil.cs)
		psi_w_stem = self.psi_wf_stem(self.vw_stem)
		psi_w_osm_stem, psi_w_turgor_stem = self.psi_components_stem(self.vw_stem)
		gw_stem = self.gwf_stem(psi_w_stem)
		psi_w_leaf = self.psi_wf_leaf(self.vw_leaf)
		psi_w_osm_leaf, psi_w_turgor_leaf = self.psi_components_leaf(self.vw_leaf)
		gw_leaf = self.gwf_leaf(psi_w_leaf)

		psi_x_val = self.psi_x(evf_val, psi_l, gp, lai, gw_leaf=gw_leaf, psi_w_leaf=psi_w_leaf)
		psi_b_val = self.psi_b(soil, soil.s, self.zr, psi_l, psi_s_arr, B, root_frac_arr, gp, lai, evf_val, gw_leaf, psi_w_leaf)
		qbx_val = self.qbx(gp, psi_x_val, psi_b_val, lai)

		# Stem and leaf storage fluxes from vwf-style storage updates.
		qwf_stem = self.qwf_stem(self.vw_stem, evf_val, gp, psi_l, lai, dt, psi_x=psi_x_val)
		qwf_leaf = self.qwf_leaf(self.vw_leaf, psi_l, lai, dt)

		energy_balance = (
			phi * lai
			- self.shf(tl, ta, lai)
			- LAMBDA_W * RHO_W * evf_val / 1_000_000.0
		)
		water_balance = evf_val - qbx_val - qwf_stem - qwf_leaf
		return (water_balance, energy_balance)

	def update(self, atm, soil, photo, root_frac_arr, B, cs_arr, dt):
		self.out = least_squares(
			self.fBal,
			(-1., 292.),
			args=(
				soil, photo, atm.phi, atm.ta, atm.qa, photo.cm,
				soil.s, self.lai, self.gp, 1., self.zr,
				root_frac_arr, B, cs_arr, dt,
			),
			bounds=([-10.0, 260.0], [0.0, 330.0]),
			method='trf', ftol=3e-16, xtol=3e-16, x_scale='jac', max_nfev=2000
		)

		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cm, self.lai, 1.)

		psi_w_stem = self.psi_wf_stem(self.vw_stem)
		psi_w_osm_stem, psi_w_turgor_stem = self.psi_components_stem(self.vw_stem)
		gw_stem = self.gwf_stem(psi_w_stem)
		psi_w_leaf = self.psi_wf_leaf(self.vw_leaf)
		psi_w_osm_leaf, psi_w_turgor_leaf = self.psi_components_leaf(self.vw_leaf)
		gw_leaf = self.gwf_leaf(psi_w_leaf)

		psi_s_arr = soil.psi_s(soil.s, soil.cs)
		psi_x_val = self.psi_x(self.ev, self.psi_l, self.gp, self.lai, gw_leaf=gw_leaf, psi_w_leaf=psi_w_leaf)
		psi_b_val = self.psi_b(soil, soil.s, self.zr, self.psi_l, psi_s_arr, B, root_frac_arr, self.gp, self.lai, self.ev, gw_leaf, psi_w_leaf)
		self.qs = self.qsf(soil, soil.s, self.zr, psi_s_arr, psi_b_val, B, root_frac_arr)
		self.qbxf = self.qbx(self.gp, psi_x_val, psi_b_val, self.lai)
		self.qxlf = self.qxl(self.gp, psi_x_val, self.psi_l, self.lai)

		self.qw_stem = self.qwf_stem(self.vw_stem, self.ev, self.gp, self.psi_l, self.lai, dt, psi_x=psi_x_val)
		self.qw_leaf = self.qwf_leaf(self.vw_leaf, self.psi_l, self.lai, dt)

		# Storage state updates
		self.vw_stem = min(max(self.vwf_stem(self.vw_stem, self.ev, self.gp, self.psi_l, self.lai, dt, psi_x=psi_x_val), 1e-12), self.VWTSTEM)
		self.vw_leaf = min(max(self.vwf_leaf(self.vw_leaf, self.psi_l, dt), 1e-12), self.VWTLEAF)
		self.w_stem = self.vw_stem / self.VWTSTEM
		self.w_leaf = self.vw_leaf / self.VWTLEAF

		uptake_rate = self.uptake(self.qs, cs_arr)
		uptake_step = uptake_rate * dt
		uptake_leaf_step = self.leaf_uptake_frac * uptake_step
		uptake_stem_step = uptake_step - uptake_leaf_step

		self.MW_uptake_stem += uptake_stem_step
		self.MW_uptake_leaf += uptake_leaf_step
		self.MW_uptake_total += uptake_step

		# Track storage concentrations for legacy mode and diagnostics.
		self.MW_stem += uptake_stem_step
		self.MW_leaf += uptake_leaf_step
		# Units: vw_stem/leaf = m^3 water/m^2 leaf, lai = m^2 leaf/m^2 ground.
		# c_stem and c_leaf computed below are mol/m3 on a ground area basis.
		self.c_stem = self.MW_stem / max(self.vw_stem * self.lai, 1e-12)
		self.c_leaf = self.MW_leaf / max(self.vw_leaf * self.lai, 1e-12)

		self._update_pi0_from_uptake()

		self.flux_balance.append(np.sum(self.qs) + self.qw_stem + self.qw_leaf - self.ev)

		# Store outputs
		self.qs_a.append(self.qs)
		gsr_vals = self.gsr(soil, soil.s, self.zr, B, root_frac_arr)
		for i, gsr_val in enumerate(gsr_vals):
			self.gsr_a[i].append(gsr_val)

		self.psi_l_a.append(self.psi_l)
		self.tl_a.append(self.tl)
		self.ev_a.append(self.ev)
		self.gp_a.append(self.gp)
		self.gsw_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))

		self.psi_x_a.append(psi_x_val)
		self.psi_b_a.append(psi_b_val)
		self.qbx_a.append(self.qbxf)
		self.qxl_a.append(self.qxlf)

		self.qw_stem_a.append(self.qw_stem)
		self.gw_stem_a.append(gw_stem)
		self.psi_w_stem_a.append(psi_w_stem)
		self.psi_w_osm_stem_a.append(psi_w_osm_stem)
		self.psi_w_turgor_stem_a.append(psi_w_turgor_stem)
		self.vw_stem_a.append(self.vw_stem)
		self.c_stem_a.append(self.c_stem)
		self.pi0_stem_a.append(self.pi0_stem)
		self.psi_uptake_stem_a.append(self.psi_uptake_stem)
		self.MW_uptake_stem_a.append(self.MW_uptake_stem)

		self.qw_leaf_a.append(self.qw_leaf)
		self.gw_leaf_a.append(gw_leaf)
		self.psi_w_leaf_a.append(psi_w_leaf)
		self.psi_w_osm_leaf_a.append(psi_w_osm_leaf)
		self.psi_w_turgor_leaf_a.append(psi_w_turgor_leaf)
		self.vw_leaf_a.append(self.vw_leaf)
		self.c_leaf_a.append(self.c_leaf)
		self.pi0_leaf_a.append(self.pi0_leaf)
		self.psi_uptake_leaf_a.append(self.psi_uptake_leaf)
		self.MW_uptake_leaf_a.append(self.MW_uptake_leaf)
		self.MW_uptake_total_a.append(self.MW_uptake_total)

	def output(self):
		return {
			'psi_l': self.psi_l_a,
			'qs': self.qs_a,
			'gp': self.gp_a,
			'gsr': self.gsr_a,
			'gsw': self.gsw_a,
			'tl': self.tl_a,
			'ev': self.ev_a,
			'ev_cum': np.cumsum(list(i * 1.8 for i in self.ev_a)),
			'psi_x': self.psi_x_a,
			'psi_b': self.psi_b_a,
			'qbx': self.qbx_a,
			'qxl': self.qxl_a,
			'flux_balance': self.flux_balance,
			# Stem storage outputs (subscripted)
			'qw_stem': self.qw_stem_a,
			'gw_stem': self.gw_stem_a,
			'psi_w_stem': self.psi_w_stem_a,
			'psi_w_osm_stem': self.psi_w_osm_stem_a,
			'psi_w_turgor_stem': self.psi_w_turgor_stem_a,
			'vw_stem': self.vw_stem_a,
			'c_stem': self.c_stem_a,
			'pi0_stem': self.pi0_stem_a,
			'psi_uptake_stem': self.psi_uptake_stem_a,
			'psi_uptake_stem_cum': self.psi_uptake_stem_a,
			'MW_uptake_stem': self.MW_uptake_stem_a,
			'w_stem': [vw / self.VWTSTEM for vw in self.vw_stem_a],
			# Leaf storage outputs
			'qw_leaf': self.qw_leaf_a,
			'gw_leaf': self.gw_leaf_a,
			'psi_w_leaf': self.psi_w_leaf_a,
			'psi_w_osm_leaf': self.psi_w_osm_leaf_a,
			'psi_w_turgor_leaf': self.psi_w_turgor_leaf_a,
			'vw_leaf': self.vw_leaf_a,
			'c_leaf': self.c_leaf_a,
			'pi0_leaf': self.pi0_leaf_a,
			'psi_uptake_leaf': self.psi_uptake_leaf_a,
			'psi_uptake_leaf_cum': self.psi_uptake_leaf_a,
			'MW_uptake_leaf': self.MW_uptake_leaf_a,
			'MW_uptake_total': self.MW_uptake_total_a,
			'leaf_uptake_frac': self.leaf_uptake_frac,
			'w_leaf': [vw / self.VWTLEAF for vw in self.vw_leaf_a],
			# Compatibility aliases
			'qw': self.qw_stem_a,
			'psi_w': self.psi_w_stem_a,
			'psi_w_osm': self.psi_w_osm_stem_a,
			'psi_w_turgor': self.psi_w_turgor_stem_a,
			'gwf': self.gw_stem_a,
			'vw': self.vw_stem_a,
		}


class HalophyteNoStorageMultiComp(Hydro):
	"""Multi-compartment, no-storage hydraulics compatible with SimulationMultiComp."""

	def __init__(self, species, atm, soil, photo, vwi_stem, c_stem=None, s_arr=None, root_frac_arr=None, B=1.0,
				 cs_arr=None, wr_stem=None, wft_stem=None, pi0_stem=None, eta_stem=None, mcap_stem=None, dt=1800,
				 salt_uptake=False, psi_wf_mode='bartlett'):
		Hydro.__init__(self, species)
		# c_stem retained for interface compatibility; when provided it should be mol/m3 on a ground area basis.

		if s_arr is None:
			s_arr = soil.s
		if root_frac_arr is None:
			root_frac_arr = np.ones(len(s_arr)) / float(len(s_arr))
		if cs_arr is None and hasattr(soil, 'cs'):
			cs_arr = soil.cs

		# Keep Halophyte-like constructor compatibility while using a no-storage formulation.
		self.dt = dt
		self.salt_uptake = salt_uptake
		self.psi_wf_mode = psi_wf_mode
		self.root_frac_arr = root_frac_arr
		self.B = B

		# Time-series outputs
		self.qs_a = []
		self.psi_l_a = []
		self.tl_a = []
		self.ev_a = []
		self.gp_a = []
		self.gsw_a = []
		self.qbx_a = []
		self.psi_b_a = []

		self.out = least_squares(
			self.fBal,
			(-1., 292.),
			args=(
				soil, photo, atm.phi, atm.ta, atm.qa, self._photo_ci(photo),
				s_arr, self.lai, self.gp, 1., self.zr, root_frac_arr, B, cs_arr
			),
			bounds=([-10.0, 260.0], [0.0, 330.0]),
			method='trf', ftol=3e-16, xtol=3e-16, x_scale='jac', max_nfev=2000
		)
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, self._photo_ci(photo), self.lai, 1.)

		psi_s_arr = self._psi_s_arr(soil, s_arr, cs_arr)
		self.psi_b_val = self.psi_b(soil, s_arr, self.zr, self.psi_l, psi_s_arr, B, root_frac_arr, self.gp, self.lai)
		self.qs = self.qsf(soil, s_arr, self.zr, psi_s_arr, self.psi_b_val, B, root_frac_arr)
		self.qbx_val = self.qbx(self.gp, self.psi_b_val, self.psi_l, self.lai)

	def _photo_ci(self, photo):
		"""Support photo models that expose either cm or cx."""
		if hasattr(photo, 'cm'):
			return photo.cm
		return photo.cx

	def _psi_s_arr(self, soil, s_arr, cs_arr=None):
		"""Compute soil water potential for all compartments for salty and non-salty soils."""
		if cs_arr is not None:
			return soil.psi_s(s_arr, cs_arr)
		try:
			return soil.psi_s(s_arr)
		except TypeError:
			if hasattr(soil, 'cs'):
				return soil.psi_s(s_arr, soil.cs)
			raise

	def update(self, atm, soil, photo, root_frac_arr, B, cs_arr, dt):
		if root_frac_arr is None:
			root_frac_arr = self.root_frac_arr
		if B is None:
			B = self.B
		if cs_arr is None and hasattr(soil, 'cs'):
			cs_arr = soil.cs
		self.root_frac_arr = root_frac_arr
		self.B = B

		self.out = least_squares(
			self.fBal,
			(-1., 292.),
			args=(
				soil, photo, atm.phi, atm.ta, atm.qa, self._photo_ci(photo),
				soil.s, self.lai, self.gp, 1., self.zr, root_frac_arr, B, cs_arr
			),
			bounds=([-10.0, 260.0], [0.0, 330.0]),
			method='trf', ftol=3e-16, xtol=3e-16, x_scale='jac', max_nfev=2000
		)
		self.psi_l = self.out.x[0]
		self.tl = self.out.x[1]
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, self._photo_ci(photo), self.lai, 1.)

		psi_s_arr = self._psi_s_arr(soil, soil.s, cs_arr)
		self.psi_b_val = self.psi_b(soil, soil.s, self.zr, self.psi_l, psi_s_arr, B, root_frac_arr, self.gp, self.lai)
		self.qs = self.qsf(soil, soil.s, self.zr, psi_s_arr, self.psi_b_val, B, root_frac_arr)
		self.qbx_val = self.qbx(self.gp, self.psi_b_val, self.psi_l, self.lai)

		# Store outputs
		self.qs_a.append(self.qs)
		self.psi_l_a.append(self.psi_l)
		self.tl_a.append(self.tl)
		self.ev_a.append(self.ev)
		self.gp_a.append(self.gp)
		self.gsw_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, self._photo_ci(photo), 1.))
		self.qbx_a.append(self.qbx_val)
		self.psi_b_a.append(self.psi_b_val)

	def output(self):
		return {
			'psi_l': self.psi_l_a,
			'qs': self.qs_a,
			'gp': self.gp_a,
			'gsw': self.gsw_a,
			'tl': self.tl_a,
			'ev': self.ev_a,
			'ev_cum': np.cumsum(list(i * 1.8 for i in self.ev_a)),
			'psi_b': self.psi_b_a,
			'qbx': self.qbx_a,
		}

	def a(self, soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr):
		"""Weighted sum of soil-root conductance times soil water potential for all compartments."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return np.sum(gsr_vals * np.array(psi_s_arr))

	def d(self, soil, s_arr, zr, psi_l, B, root_frac_arr):
		"""Sum of soil-root conductances for all compartments."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return np.sum(gsr_vals)

	def psi_b(self, soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr, gp, lai):
		"""Basal root node potential without storage compartment."""
		a_val = self.a(soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr)
		d_val = self.d(soil, s_arr, zr, psi_l, B, root_frac_arr)
		gp_term = gp * lai
		denom = d_val + gp_term
		if abs(denom) < 1e-12:
			denom = 1e-12 if denom >= 0 else -1e-12
		return (a_val + gp_term * psi_l) / denom

	def qsf(self, soil, s_arr, zr, psi_s_arr, psi_b, B, root_frac_arr):
		"""Soil water flux for multiple compartments (array output)."""
		gsr_vals = self.gsr(soil, s_arr, zr, B, root_frac_arr)
		return gsr_vals * (np.array(psi_s_arr) - psi_b)

	def qbx(self, gp, psi_b, psi_l, lai):
		"""Basal-to-xylem flow for no-storage configuration."""
		return gp * lai * (psi_b - psi_l)

	def gsr(self, soil, s_arr, zr, B, root_frac_arr):
		"""Soil-Root conductance for multiple compartments (array output, B is constant)."""
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

	def fBal(self, params, soil, photo, phi, ta, qa, c1, s_arr, lai, gp, ared, zr, root_frac_arr, B, cs_arr):
		"""Energy and water balance equations for multi-compartment Halophyte without storage."""
		psi_l, tl = params
		evf_val = self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)
		psi_s_arr = self._psi_s_arr(soil, s_arr, cs_arr)
		psi_b_val = self.psi_b(soil, s_arr, zr, psi_l, psi_s_arr, B, root_frac_arr, gp, lai)
		qbx_val = self.qbx(gp, psi_b_val, psi_l, lai)

		energy_balance = (
			phi * lai
			- self.shf(tl, ta, lai)
			- LAMBDA_W * RHO_W * evf_val / 1_000_000.0
		)
		water_balance = evf_val - qbx_val
		return (water_balance, energy_balance)
