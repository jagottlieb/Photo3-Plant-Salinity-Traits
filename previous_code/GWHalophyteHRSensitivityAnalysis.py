# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 21:07:00 2022

@author: joshg
"""

from scipy.optimize import *
#from sympy import *
import numpy as np
import pandas as pd
from math import exp, pi, sqrt, log
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('TkAgg')
from dics import *
from functions import *
from defs import *
from GWHalophyteHR import *
import openpyxl
import time

s_list = list(np.linspace(0.09,0.21, 40))
cs_list = list(np.linspace(0, 260, 40))

HR_LIST_m = [] #List for 5 day cumulative HR results for each iteration of cw
HR_LIST_o = [] #List for 5 day cumulative HR results for each iteration of cw
psi_s_mat_list = [] #List of accompanying psi_s matric values
m_gradient_list = [] #List of initial matric potential gradient between GW and Soil 
psi_s_o_list = [] #List of accompanying psi_s osmotic values
o_gradient_list = [] #List of initial osmotic potential gradient between GW and Soil

for s in s_list:
	duration = 5
	weatherFile = "sample_data\Almeria2017July.xlsx"
	resultsFile = 'sample_output\SolverTrials\TrialOctober10.csv' # default value for the location where results are saved
	
	timestepM = 30 # Model change in time at each step (min)
	timestepD = 30 # timestep of input data 
	dt = timestepM*60. # no. of seconds in timestep, used to advance differential equations
	df = pd.read_excel(weatherFile, engine='openpyxl')
	tempC = df['Temperature'] # reads in C
	taInp = tempC + 273. # convert to K
	rh = df['Relative Humidity'] # extracts rh column (%)
	psat = A_SAT*np.exp((B_SAT*(tempC))/(C_SAT + tempC)) # saturated vapor pressure in Pa
	qaInp = 0.622*rh/100.*psat/P_ATM # needs to be in kg/kg
	qaInp = list(qaInp.values)
	taInp = list(taInp.values)
	phiInp = list(df['GHI'].values)  # extracts global solar radiation column from Excel Worksheet in W/m^2
	
	# enter values manually
	sinit = s #0.12 #Initial soil moisture
	vwi = 0.86 #0.81 #Initial plant storage volume
	gw_z_init = 3 #Initial depth to the groundwater table (m)
	cs_init = 200 #200 #220 #5 #100 #initial salt concentration in soil, (Mol/m^3)
	cgw_init = 300 #300 #345 (from TDS conversion of Armas data) ##initial salt concentration in groundwater (Mol/m^3)
	cw_init = 300 #50 #initial salt concentration in plant storage (Mol/m^3)
	species = LentiscPistacia() # Mastic Tree, Lentisc pistacia
	
	#Classes for simulations
	atmosphere = Atmosphere(phiInp[0], taInp[0], qaInp[0])
	soil = SaltySoilGW(Sand(), DrydownSoil(), species.ZR, sinit, cs_init, cgw_init, gw_z_init)
	photo = C3(species, atmosphere)
	hydro = HalophyteGW(species, atmosphere, soil, photo, vwi, cw_init)
	plant = Simulation(species, atmosphere, soil, photo, hydro)
	
	
	for i in range(steps(duration, int(timestepM))):
		plant.update(dt, phiInp[i], taInp[i], qaInp[i])
	results = plant.output()
	
	HR_LIST_m.append(results['hr_cum'][-1])
	psi_s_mat_list.append(results['psi_s_mat'][0])
	m_gradient_list.append(results['psi_gw_mat'][0]-results['psi_s_mat'][0])
	
for cs in cs_list:
	duration = 5
	weatherFile = "sample_data\Almeria2017July.xlsx"
	resultsFile = 'sample_output\SolverTrials\TrialOctober10.csv' # default value for the location where results are saved
	
	timestepM = 30 # Model change in time at each step (min)
	timestepD = 30 # timestep of input data 
	dt = timestepM*60. # no. of seconds in timestep, used to advance differential equations
	df = pd.read_excel(weatherFile, engine='openpyxl')
	tempC = df['Temperature'] # reads in C
	taInp = tempC + 273. # convert to K
	rh = df['Relative Humidity'] # extracts rh column (%)
	psat = A_SAT*np.exp((B_SAT*(tempC))/(C_SAT + tempC)) # saturated vapor pressure in Pa
	qaInp = 0.622*rh/100.*psat/P_ATM # needs to be in kg/kg
	qaInp = list(qaInp.values)
	taInp = list(taInp.values)
	phiInp = list(df['GHI'].values)  # extracts global solar radiation column from Excel Worksheet in W/m^2
	
	# enter values manually
	sinit = 0.12 #0.12 #Initial soil moisture
	vwi = 0.86 #0.81 #Initial plant storage volume
	gw_z_init = 3 #Initial depth to the groundwater table (m)
	cs_init = cs #200 #220 #5 #100 #initial salt concentration in soil, (Mol/m^3)
	cgw_init = 345 #300 #345 (from TDS conversion of Armas data) ##initial salt concentration in groundwater (Mol/m^3)
	cw_init = 270 #50 #initial salt concentration in plant storage (Mol/m^3)
	species = LentiscPistacia() # Mastic Tree, Lentisc pistacia
	
	#Classes for simulations
	atmosphere = Atmosphere(phiInp[0], taInp[0], qaInp[0])
	soil = SaltySoilGW(Sand(), DrydownSoil(), species.ZR, sinit, cs_init, cgw_init, gw_z_init)
	photo = C3(species, atmosphere)
	hydro = HalophyteGW(species, atmosphere, soil, photo, vwi, cw_init)
	plant = Simulation(species, atmosphere, soil, photo, hydro)
	
	for i in range(steps(duration, int(timestepM))):
		plant.update(dt, phiInp[i], taInp[i], qaInp[i])
	results = plant.output()
	
	HR_LIST_o.append(results['hr_cum'][-1])
	psi_s_o_list.append(results['psi_s_osm'][0])
	o_gradient_list.append(results['psi_gw_osm'][0]-results['psi_s_osm'][0])
	print(HR_LIST_o)


hr_mag = plt.figure()
hr_mag.suptitle('5 day HR vs Soil Potential Gradients')
plt.xlabel('Potential Gradient (MPa)')
plt.ylabel('HR (mm)')
#plt.xlim([psi_gradient_list[0],psi_s_mat_list[-1]])
#plt.ylim([0,max(HR_LIST_m)+1])
plt.plot(m_gradient_list, HR_LIST_m, 'k-', label = 'Matric')
plt.plot(o_gradient_list, HR_LIST_o, 'r-', label ='Osmotic')
#plt.minorticks_on()
plt.legend('Matric', 'Osmotic')
plt.grid(b ='true', which = 'both',axis = 'both')
hr_mag.show()

# 	data = pd.DataFrame.from_dict(results)
# 	datacsv = pd.DataFrame.from_dict(results)
# 	datacsv.to_csv(resultsFile)
	#data.to_pickle(resultsFile)

	
	
	#Plot results
# 	startDay = 0
# 	endDay = duration
# 	dispDuration = endDay-startDay
# 	daySteps = 60//timestepM*24
# 	timevec = np.linspace(0,duration,duration*daySteps)
# 	timevecHr = np.linspace(0,duration*24,duration*daySteps)
# 	
# 	
# 	anp = plt.figure()
# 	gph_title = 'Soil moisture: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init)
# 	plt.title(gph_title)
# 	plt.xlabel("time (days)")
# 	plt.ylabel("s (-)")
# 	plt.plot(timevec, results['s'][daySteps*startDay:daySteps*endDay])
# 	for i in range(0, endDay):
# 	 	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
# 	plt.xlim(0, dispDuration)
# 	plt.ylim(0, 0.3)
# 	plt.xticks([0.,6.,12.,18.,24.,30.])
# 	plt.grid(true, axis='y')
# 	plt.legend()
# 	anp.show()
# 	
# 	fluxes = plt.figure()
# 	fluxes.suptitle('Fluxes: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
# 	plt.xlabel("t (days)")
# 	plt.ylabel("Flux (um/s)")
# 	plt.ylim([-0.01,0.1])
# 	plt.xlim([0,endDay])
# 	print(timevec)
# 	plt.plot(timevec, results['qs'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
# 	#plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['qs'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
# 	plt.plot(timevec, results['qgw'][daySteps*startDay:daySteps*endDay], 'k-', label = 'GW')
# 	plt.plot(timevec, results['ev'][daySteps*startDay:daySteps*endDay], 'b:', label = 'ET')
# 	plt.plot(timevec, results['qw'][daySteps*startDay:daySteps*endDay], 'g:', label = 'Storage')
# 	plt.legend(['Soil','GW','EV', 'Storage'], loc='upper right')
# 	plt.minorticks_on()
# 	plt.grid(b ='true', which = 'both', axis = 'both')
# 	
# 	potents = plt.figure()
# 	potents.suptitle('Potentials: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
# 	plt.xlabel('time (days)')
# 	plt.ylabel('Potential (MPa)')
# 	plt.xlim([0, endDay])
# 	plt.ylim([-3,0.5])
# 	plt.plot(timevec, results['psi_s'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
# 	plt.plot(timevec, results['psi_b'][daySteps*startDay:daySteps*endDay], 'k-', label = 'Base')
# 	plt.plot(timevec, results['psi_x'][daySteps*startDay:daySteps*endDay], 'b:', label = 'Node')
# 	plt.plot(timevec, results['psi_l'][daySteps*startDay:daySteps*endDay], 'g-', label = 'Leaf')
# 	plt.plot(timevec, results['psi_gw'][daySteps*startDay:daySteps*endDay], 'r:', label = 'GW')
# 	plt.plot(timevec, results['psi_w'][daySteps*startDay:daySteps*endDay], 'c:', label = 'PWS')
# 	plt.legend(['Soil','Base','Node', 'Leaf', 'GW', 'PWS'])
# 	plt.minorticks_on()
# 	plt.grid(b ='true', which = 'both', axis = 'both')
# 	
# 	##### Potential Plots of Plant, Soil, Subsoil
# 	total_potent = plt.figure()
# 	total_potent.suptitle('Total Potential: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
# 	plt.xlabel('time (days)')
# 	plt.ylabel('Potential (MPa)')
# 	plt.xlim([0, endDay])
# 	plt.ylim([-3,0.5])
# 	plt.plot(timevec, results['psi_s'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
# 	plt.plot(timevec, results['psi_gw'][daySteps*startDay:daySteps*endDay], 'r:', label = 'GW')
# 	plt.plot(timevec, results['psi_b'][daySteps*startDay:daySteps*endDay], 'c:', label = 'Base')
# 	plt.legend(['Soil','GW', 'Base'])
# 	plt.minorticks_on()
# 	plt.grid(b ='true', which = 'both',axis = 'both')
# 	
# 	osm_potent = plt.figure()
# 	osm_potent.suptitle('Osmotic Potential: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
# 	plt.xlabel('time (days)')
# 	plt.ylabel('Potential (MPa)')
# 	plt.xlim([0, endDay])
# 	plt.ylim([-3,0.5])
# 	plt.plot(timevec, results['psi_s_osm'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
# 	plt.plot(timevec, results['psi_gw_osm'][daySteps*startDay:daySteps*endDay], 'r:', label = 'GW')
# 	plt.plot(timevec, results['psi_w_osm'][daySteps*startDay:daySteps*endDay], 'c:', label = 'PWS')
# 	plt.legend(['Soil','GW', 'PWS'])
# 	plt.minorticks_on()
# 	plt.grid(b ='true', which = 'both',axis = 'both')
# 	
# 	mat_potent = plt.figure()
# 	mat_potent.suptitle('Matric/Turgor Potential: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
# 	plt.xlabel('time (days)')
# 	plt.ylabel('Potential (MPa)')
# 	plt.xlim([0, endDay])
# 	plt.ylim([-3,1.3])
# 	plt.plot(timevec, results['psi_s_mat'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
# 	plt.plot(timevec, results['psi_gw_mat'][daySteps*startDay:daySteps*endDay], 'r:', label = 'GW')
# 	plt.plot(timevec, results['psi_w_turgor'][daySteps*startDay:daySteps*endDay], 'c:', label = 'PWS')
# 	plt.legend(['Soil','GW', 'PWS'])
# 	plt.minorticks_on()
# 	plt.grid(b ='true', which = 'both',axis = 'both')
# 	
# 	hr_cum = plt.figure()
# 	hr_cum.suptitle('Cumulative HR: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
# 	plt.xlabel('time (days)')
# 	plt.ylabel('Cumulative HR (mm)')
# 	plt.xlim([0, endDay])
# 	#plt.ylim([-3,0.5])
# 	plt.plot(timevec, results['hr_cum'][daySteps*startDay:daySteps*endDay], 'm-')
# 	plt.minorticks_on()
# 	plt.grid(b ='true', which = 'both',axis = 'both')
# 	
# 	hr_matric = plt.figure()
# 	hr_matric.suptitle('Cumulative HR vs Soil Matric Potential: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
# 	plt.xlabel('Matric Potential (-MPa)')
# 	plt.ylabel('Cumulative HR (mm)')
# 	#plt.xlim([0, endDay])
# 	#plt.ylim([-3,0.5])
# 	plt.plot([i*(-1) for i in results['psi_s_mat'][daySteps*startDay:daySteps*endDay]], results['hr_cum'][daySteps*startDay:daySteps*endDay], 'k-')
# 	plt.minorticks_on()
# 	plt.grid(b ='true', which = 'both',axis = 'both')

# 	startDay = 0
# 	endDay = duration
# 	dispDuration = endDay-startDay
# 	daySteps = 60//timestepM*24
# 	timevec = np.linspace(0,duration,duration*daySteps)
# 	timevecHr = np.linspace(0,duration*24,duration*daySteps)