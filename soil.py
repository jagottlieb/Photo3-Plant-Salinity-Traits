from scipy.optimize import *
from dics import *
from functions import *
import numpy as np
import time

class DrydownSoil(object):
	def __init__(self):
		pass
	def snew(self, soil, dt, zr, qs):
		new_s = (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s
		# Clamp to physical range to prevent negative moisture
		return float(np.clip(new_s, 1e-9, 1.0))

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
		self.stype = stype
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
		# Use the texture's own retention curve (Van Genuchten for Berger, Clapp-Hornberger otherwise)
		if getattr(self.stype, 'VAN_GENUCHTEN', False):
			matric = self.stype.psi_s(s)
		else:
			matric = self.PSI_SS*(s**-self.B)
		return matric - self.E*self.cs*R*self.IV*self.TS*10.**(-6.)

# --- Multi-compartment classes ---

class SoilMultiple(object):
	EVMAX = 3
	SY = 0.2
	def __init__(self, stype, dynamics, zr, s, cs=None):
		self.s = np.atleast_1d(s)
		self.ZR = np.atleast_1d(zr)
		# stype can be array or single
		if isinstance(stype, (list, tuple, np.ndarray)):
			if len(stype) == 1 and len(self.s) > 1:
				self.stype = [stype[0] for _ in self.s]
			else:
				self.stype = list(stype)
		else:
			self.stype = [stype for _ in self.s]
		# Dynamics can be a single object or list
		if isinstance(dynamics, (list, tuple, np.ndarray)):
			if len(dynamics) == 1 and len(self.s) > 1:
				self.dynamics = [dynamics[0] for _ in self.s]
			else:
				self.dynamics = list(dynamics)
		else:
			self.dynamics = [dynamics for _ in self.s]
		# Broadcast ZR if needed
		if len(self.ZR) == 1 and len(self.s) > 1:
			self.ZR = np.array([self.ZR[0] for _ in self.s])
		# cs can be array or single or None
		if cs is not None:
			self.cs = np.atleast_1d(cs)
			if len(self.cs) == 1 and len(self.s) > 1:
				self.cs = np.array([self.cs[0] for _ in self.s])
		else:
			self.cs = None
		self.PSI_SS = self.stype[0].PSI_SS
		self.B = self.stype[0].B
		self.KS = self.stype[0].KS
		self.N = self.stype[0].N
		self.SH = self.stype[0].SH
		self.s_a = [[] for _ in self.s]
		self.psi_s_a = [[] for _ in self.s]
		self.cs_a = [[] for _ in self.s] if self.cs is not None else None
		self.rain_amt = 0
		self.sm_inp = self.s
	def update(self, dt, zr, qs):
		qs = np.atleast_1d(qs)
		zr = np.atleast_1d(zr)
		# Store the array temporarily
		s_array = self.s.copy()
		for i in range(len(self.s)):
			# Temporarily set soil.s to scalar for this compartment
			self.s = s_array[i]
			# Call snew with scalar s value
			new_s = self.dynamics[i].snew(self, dt, zr[i], qs[i])
			new_s = float(np.clip(new_s, 1e-9, 1.0))
			# Restore array and update compartment
			self.s = s_array
			self.s[i] = new_s
			s_array = self.s.copy()
			
			self.s_a[i].append(self.s[i])
			self.psi_s_a[i].append(self.psi_s(self.s[i], cs=self.cs[i] if self.cs is not None else None, i=i))
			if self.cs is not None:
				# SaltySoil dynamics: update cs
				if hasattr(self.dynamics[i], 'update_cs'):
					self.cs[i] = self.dynamics[i].update_cs(self, i)
				self.cs_a[i].append(self.cs[i])
	def output(self):
		out = {}
		out['s'] = self.s_a
		out['psi_s'] = self.psi_s_a
		if self.cs is not None:
			out['cs'] = self.cs_a
		return out
	def leak(self, s):
		"""Calculate soil leakage. Handles both scalar and array inputs."""
		s = np.atleast_1d(s)
		result = np.zeros_like(s, dtype=float)
		# Only compute on physically valid moisture values to avoid invalid powers
		mask = (s > 0.0) & (s < 1.0)
		if len(s) == len(self.stype):
			for j in range(len(s)):
				if mask[j]:
					result[j] = self._texture_leak(self.stype[j], s[j])
		else:
			result[mask] = [self._texture_leak(self.stype[0], sv) for sv in s[mask]]
		return result if len(result) > 1 else result[0]
	# psi_s (MPa)
	def psi_s(self, s, cs=None, i=0):
		"""Calculate soil water potential. 
		Args:
			s: soil moisture (scalar or array)
			cs: salt concentration (optional, scalar or array - if provided overrides i)
			i: compartment index (used if cs not provided and s is scalar)
		Returns:
			Soil water potential (scalar or array matching input shape)
		"""
		s = np.atleast_1d(s)
		# Prevent invalid powers for dry/negative values
		s_safe = np.clip(s, 1e-9, None)

		# Matric potential: each compartment uses its own texture's retention curve
		# (Clapp-Hornberger by default, Van Genuchten for textures like Berger that define one)
		if len(s_safe) == len(self.stype):
			matric = np.array([self._texture_psi_s(self.stype[j], s_safe[j]) for j in range(len(s_safe))])
		else:
			matric = np.array([self._texture_psi_s(self.stype[i], sv) for sv in s_safe])

		# If SaltySoil, include osmotic term
		if self.cs is not None:
			if cs is not None:
				# cs provided explicitly (can be scalar or array)
				cs_val = np.atleast_1d(cs)
			elif len(s) == len(self.cs):
				# s is array matching number of compartments, use all cs values
				cs_val = self.cs
			else:
				# s is scalar or single value, use indexed cs
				cs_val = self.cs[i]
			
			psi_s_val = matric - cs_val*R*SaltySoil.IV*SaltySoil.TS*10.**(-6.)
		else:
			psi_s_val = matric
		
		# Return scalar if input was scalar, array otherwise
		return psi_s_val if len(psi_s_val) > 1 or isinstance(cs, np.ndarray) else psi_s_val[0]
	def _texture_psi_s(self, stype, s):
		"""Matric potential (MPa) for one texture: Van Genuchten if defined, else Clapp-Hornberger."""
		if getattr(stype, 'VAN_GENUCHTEN', False):
			return stype.psi_s(s)
		return stype.PSI_SS*(s**-stype.B)
	def _texture_leak(self, stype, s):
		"""Leakage/hydraulic conductivity for one texture: Mualem-Van Genuchten if defined, else Clapp-Hornberger."""
		if getattr(stype, 'VAN_GENUCHTEN', False):
			return stype.leak(s)
		return .11574*stype.KS*s**(2.*stype.B + 3.)
	def evap(self, s):
		"""Calculate soil evaporation. Handles both scalar and array inputs."""
		return np.where(s > self.SH, self.EVMAX*(s - self.SH)/(1. - self.SH), 0.)

class SaltySoilMultiple(SoilMultiple):
	TS = 293.
	IV = 2.
	E = 0.92
	SY = 0.5
	def __init__(self, stype, dynamics, zr, s, cs):
		super().__init__(stype, dynamics, zr, s)
		self.cs = np.atleast_1d(cs)
		if len(self.cs) == 1 and len(self.s) > 1:
			self.cs = np.array([self.cs[0] for _ in self.s])
		# Salt mass per compartment, mol/m2
		self.MS = self.cs * self.ZR * self.N * self.s
		self.cs_a = [[] for _ in self.s]
		self.MS_a = [[] for _ in self.s]
	def update(self, dt, zr, qs):
		qs = np.atleast_1d(qs)
		zr = np.atleast_1d(zr)

		# Scalar-safe compartment updates (matches DrydownSoil.snew signature).
		s_array = self.s.copy()
		for i in range(len(self.s)):
			# Temporarily set soil.s to scalar for this compartment
			self.s = s_array[i]
			# Call snew with scalar s value
			new_s = self.dynamics[i].snew(self, dt, zr[i], qs[i])
			new_s = float(np.clip(new_s, 1e-9, 1.0))
			# Restore array and update compartment
			self.s = s_array
			self.s[i] = new_s
			s_array = self.s.copy()

		# Match SaltySoil mass-balance update for salt concentration, per compartment.
		self.cs = self.MS/(self.s*self.N*self.ZR)
		psi_s_vals = np.atleast_1d(self.psi_s(self.s, cs=self.cs))

		for i in range(len(self.s)):
			self.s_a[i].append(self.s[i])
			self.psi_s_a[i].append(psi_s_vals[i])
			self.cs_a[i].append(self.cs[i])
			self.MS_a[i].append(self.MS[i])

	def output(self):
		out = super().output()
		out['MS'] = self.MS_a
		return out

	def psi_s(self, s, cs=None, i=0):
		s = np.atleast_1d(s)
		s_safe = np.clip(s, 1e-9, None)
		if len(s_safe) == len(self.stype):
			matric = np.array([self._texture_psi_s(self.stype[j], s_safe[j]) for j in range(len(s_safe))])
		else:
			matric = np.array([self._texture_psi_s(self.stype[i], sv) for sv in s_safe])
		if cs is not None:
			cs_val = np.atleast_1d(cs)
		elif len(s_safe) == len(self.cs):
			cs_val = self.cs
		else:
			cs_val = self.cs[i]
		psi_s_val = matric - cs_val*R*self.IV*self.TS*10.**(-6.)
		return psi_s_val if len(psi_s_val) > 1 or isinstance(cs, np.ndarray) else psi_s_val[0]

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

class Berger(object):
	"""Berger soil texture: Van Genuchten retention curve, Mualem unsaturated hydraulic conductivity."""
	VAN_GENUCHTEN = True
	TH_R = 0.155
	TH_S = 0.807
	ALPHA = 0.0502   # 1/cm
	VG_N = 1.768
	VG_M = 1 - 1/VG_N
	L = 0.5   # pore-connectivity/tortuosity parameter
	N = 0.8   # measured porosity; relative saturation s = theta/N
	# Unused (Van Genuchten/Mualem used instead); kept so generic soil-type lookups don't fail
	PSI_SS = None
	B = None
	KS = 1.
	SH = .15
	def __init__(self):
		pass
	def psi_s(self, s):
		"""Van Genuchten matric potential (MPa) from relative soil moisture s = theta/N."""
		s = np.asarray(s, dtype=float)
		theta = np.clip(s, 1e-9, None) * self.N
		Se = np.clip((theta - self.TH_R) / (self.TH_S - self.TH_R), 1e-6, 1. - 1e-6)
		h_cm = (1. / self.ALPHA) * (Se**(-1. / self.VG_M) - 1.)**(1. / self.VG_N)
		h_kpa = h_cm * 0.0980665  # cm H2O -> kPa
		return -h_kpa * 1e-3  # kPa -> MPa, negative = suction
	def leak(self, s):
		"""Mualem-Van Genuchten unsaturated hydraulic conductivity: K(Se) = Ks*Se^L*[1-(1-Se^(1/m))^m]^2."""
		s = np.asarray(s, dtype=float)
		theta = np.clip(s, 1e-9, None) * self.N
		Se = np.clip((theta - self.TH_R) / (self.TH_S - self.TH_R), 1e-6, 1. - 1e-6)
		k_rel = Se**self.L * (1. - (1. - Se**(1. / self.VG_M))**self.VG_M)**2
		return .11574*self.KS*k_rel
