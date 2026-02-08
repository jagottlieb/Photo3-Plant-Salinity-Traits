# Halophyte hydraulic classes extracted from your GWHalophyteHR variant
# Minimal standalone module to provide SoilGW, SaltySoilGW, Hydro, HalophyteGW
from scipy.optimize import *
from dics import *
from functions import *
import numpy as np

class SoilGW(object):
	EVMAX = 3
	SY = 0.2
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
		self.rain_amt = 0
		self.sm_inp = s
	def update(self, dt, zr, qs):
		self.s = self.dynamics.snew(self, dt, zr, qs)
		self.s_a.append(self.s)
		self.psi_s_a.append(self.psi_s(self.s))
	def output(self):
		return {'s': self.s_a, 'psi_s': self.psi_s_a}
	def leak(self, s):
		try:
			ans = .11574*self.KS*s**(2.*self.B + 3.)
		except OverflowError:
			ans = 0.
		return ans
	def psi_s(self, s):
		return self.PSI_SS*(s**-self.B)
	def evap(self, s):
		if s > self.SH:
			return self.EVMAX*(s - self.SH)/(1. - self.SH)
		else:
			return 0.

class SaltySoilGW(SoilGW):
	TS = 293.
	GWTS = 293
	IV = 2.
	E = 0.95
	SY = 0.5
	def __init__(self, stype, dynamics, zr, s, cs, cgw, gw_z):
		SoilGW.__init__(self, stype, dynamics, zr, s)
		self.s = s
		self.cs = cs
		self.MS = cs*self.ZR*self.N*s
		self.cs_a = []
		self.gw_z = gw_z
		self.gw_cs = cgw
		self.S_y = self.SY
		self.psi_s_a = []
		self.gw_z_a = []
		self.gw_cs_a = []
		self.psi_s_mat_a = []
		self.psi_s_osm_a = []
	def update(self, dt, zr, qs, qgw):
		self.gw_z_a.append(self.gw_z)
		self.psi_s_a.append(self.psi_s(self.s))
		self.s_a.append(self.s)
		self.s = (dt/(self.N*zr*10.**6)*(-qs - (self.evap(self.s)*1000.)/(24.*60*60)- self.leak(self.s))) + self.s
		self.cs = self.MS/(self.s*self.N*self.ZR)
		self.gw_z = self.gw_z - dt*60*self.leak(self.s)/10**6 + qgw*60*dt/10**6
		self.gw_cs_a.append(self.gw_cs)
		self.cs_a.append(self.cs)
		self.psi_gw_a.append(self.psi_gw(self.gw_cs))
		self.psi_s_mat_a.append(self.psi_s_mat(self.s))
		self.psi_s_osm_a.append(self.psi_s_osm(self.s))
	def output(self):
		return {'s': self.s_a, 'cs': self.cs_a,'gw_z': self.gw_z_a, 'psi_s': self.psi_s_a}
	def psi_s(self, s):
		return self.PSI_SS*(s**-self.B) - self.E*self.cs*R*self.IV*self.TS*10.**(-6.)
	def psi_s_mat(self, s):
		return self.PSI_SS*(s**-self.B)
	def psi_s_osm(self, s):
		return -self.E*self.cs*R*self.IV*self.TS*10.**(-6.)
	def psi_gw(self, gw_cs):
		return self.PSI_SS - self.E*gw_cs*R*self.IV*self.GWTS*10.**(-6.)

class Hydro(object):
	A_ROOT_l = 8.
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
		self.qgw_a = []
		self.qs_a = []
		self.qbx_a = []
	def gsw(self, photo, phi, ta, psi_l, qa, tl, ci, ared):
		return photo.gsc(phi, ta, psi_l, qa, tl, ci, ared)*1.6 + (self.GCUT*P_ATM/(1000.*R*ta))

class HalophyteGW(Hydro):
	F_CAP = 0.5
	E = 0.95
	TS = 293.
	GWTS = 293
	IV = 2.
	def __init__(self, species, atm, soil, photo, vwi, cw):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.vw = vwi*self.VWT
		self.MW = cw*self.vw
		self.CAP = species.CAP
		self.cw = cw
		self.w = self.vw/self.VWT
		self.zgw = 3.5
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
	def psi_wf(self, vw, cw):
		TL = 293.
		return (vw/self.VWT - 0.028)**8 - cw*R*self.IV*TL*10.**(-6.)
	def gwf(self, psi_w):
		return self.GWMAX*(self.vw/self.VWT)**4
	def vwf(self, vw, ev, gp, psi_l, lai, cw, dt):
		psi_w = self.psi_wf(vw, cw)
		return min(vw - self.gwf(psi_w)*(psi_w - (ev*(1. - self.F_CAP))/(lai*gp) - psi_l)*dt/10.**6, self.VWT)
	def qwf(self, vw, ev, gp, psi_l, lai, cw, dt):
		return (vw - self.vwf(vw, ev, gp, psi_l, lai, cw, dt))*lai*10.**6/dt
