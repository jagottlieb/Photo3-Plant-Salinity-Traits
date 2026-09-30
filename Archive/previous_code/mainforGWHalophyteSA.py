from genericpath import exists
import os
from scipy.optimize import * #fsolve
from sympy import *
import numpy as np
import pandas as pd
from math import exp, pi, sqrt, log
import matplotlib.pyplot as plt
from dics import *
from functions import *
import importlib as importlib
import defs
importlib.reload(defs)
from defs import *
from GWHalophyteFixingStorageEqn import *
#from GWHalophyteHRStorageFE import *
#from GWHalophyteHRSensitivityAnalysis import *
import openpyxl

for x in range(0, 1):

	duration = 20
	weatherFile = r"sample_data\AgadirInterp30.xlsx"
	resultsFolder = 'sample_output\\Diagnostic_020726_GWHalophyteSSensitivityAnalysis\\'
	resultsFile = resultsFolder + 'AgadirGWHalophyte_Diagnostic_020726.csv' # default value for the location where results are saved
	
	# Create results folder if it doesn't exist
	if not os.path.exists(resultsFolder):
		os.makedirs(resultsFolder)
	
	timestepM = 30 # Model change in time at each step (min)
	timestepD = 30 # timestep of input data 
	dt = timestepM*60. # no. of seconds in timestep, used to advance differential equations
	#plant_s = [0, 0, 0, 50, 50, 50, 150, 150, 150, 250, 250, 250, 350, 350, 350]
	#gw_s = [0, 100, 200, 0, 100, 200, 0, 100, 200, 0, 100, 200, 0, 100, 200]
	#s_s = [100, 200, 300, 100, 200, 300, 100, 200, 300, 100, 200, 300, 100, 200, 300]
	
	#gw_s = [0, 100, 200, 0, 100, 200, 0, 100, 200, 0, 100, 200, 0, 100, 200]
	#s_s = [200, 300, 400, 200, 300, 400, 200, 300, 400, 200, 300, 400, 200, 300, 400]
	
	gw_s = [150,]
	s_s = [150,]
	plant_s = [150,]
	lam_s_x = [4000,]   #Soil root length densities for each simulation
	lam_gw_x = [4000,]  #Tap root length densities for each simulation


	
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
	sinit = 0.5 #Initial soil moisture
	vwi = .95 #Initial storage volume
	gw_z_init = 3 #Initial depth to the groundwater table (m)
	cs_init = s_s[x] #initial salt concentration in soil, (Mol/m^3)
	cgw_init = gw_s[x] #initial salt concentration in groundwater (Mol/m^3)
	cw_init = plant_s[x] #initial salt concentration in plant storage (Mol/m^3)
	species = Pmenz() # Douglas fir, Pseudotsuga menziesii
	atmosphere = Atmosphere(phiInp[0], taInp[0], qaInp[0])
	soil = SaltySoilGW(Loam(), DrydownSoil(), species.ZR, sinit, cs_init, cgw_init, gw_z_init)
	photo = C3(species, atmosphere)
	#hydro = HalophyteGW(species, atmosphere, soil, photo, vwi, cw_init, lam_s_x[x], lam_gw_x[x]) # When using GWHalophyteFixingStorageEqn, use this line and comment out the line below
	hydro = HalophyteGW(species, atmosphere, soil, photo, vwi, cw_init) # When using GWHalophyteOsmoregulation, use this line and comment out the line above
	
	
	
	plant = Simulation(species, atmosphere, soil, photo, hydro)
	
	
	for i in range(steps(duration, int(timestepM))):
	
		plant.update(dt, phiInp[i], taInp[i], qaInp[i])
	
	results = plant.output()
	
	data = pd.DataFrame.from_dict(results)
	datacsv =pd.DataFrame.from_dict(results)
	datacsv.to_csv(resultsFile)
	#data.to_pickle(resultsFile)
	# Save data as csv file
	#data.to_csv(resultsFile)
	
	#Plot results
	startDay = 0
	endDay = duration
	dispDuration = endDay-startDay
	daySteps = 60//timestepM*24
	timevec = np.linspace(0,duration,duration*daySteps)
	timevecHr = np.linspace(0,duration*24,duration*daySteps)
	
	anp = plt.figure()
	gph_title = ('Soil moisture: cw {} cs {} cgw {} rld_s{} rld_gw{}'.format(plant_s[x], s_s[x], gw_s[x], lam_s_x[x], lam_gw_x[x]))
	plt.title(gph_title)
	plt.xlabel("time (days)")
	plt.ylabel("s (-)")
	plt.plot(timevec[0:daySteps*dispDuration], results['s'][daySteps*startDay:daySteps*endDay])
	for i in range(0, 25):
		plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
	plt.xlim(0, dispDuration)
	plt.ylim(0.35, 0.8)
	plt.xticks([0.,6.,12.,18.,24.])
	plt.grid(visible=True, axis='y')
	plt.legend()
	
	fluxes = plt.figure()
	fluxes.suptitle('Fluxes: cw {} cs {} cgw {} rld_s{} rld_gw{}'.format(plant_s[x], s_s[x], gw_s[x], lam_s_x[x], lam_gw_x[x]))
	plt.xlabel("t (h)")
	plt.ylabel("Flux (um/s)")
	plt.ylim([-0.1,0.6])
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['qs'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['qgw'][daySteps*startDay:daySteps*endDay], 'k-', label = 'GW')
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['ev'][daySteps*startDay:daySteps*endDay], 'b:', label = 'ET')
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['qw'][daySteps*startDay:daySteps*endDay], 'g:', label = 'Storage')
	plt.legend(['Soil','GW','EV','Storage'])
	plt.minorticks_on()
	plt.grid(visible=True, which='both', axis='both')
	
	potents = plt.figure()
	potents.suptitle('Potentials: cw {} cs {} cgw {} rld_s{} rld_gw{}'.format(plant_s[x], s_s[x], gw_s[x], lam_s_x[x], lam_gw_x[x]))
	plt.xlabel('time (h)')
	plt.ylabel('Potential (MPa)')
	plt.ylim([-3,0.5])
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['psi_s'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['psi_b'][daySteps*startDay:daySteps*endDay], 'k-', label = 'Base')
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['psi_x'][daySteps*startDay:daySteps*endDay], 'b:', label = 'Node')
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['psi_l'][daySteps*startDay:daySteps*endDay], 'g-', label = 'Leaf')
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['psi_gw'][daySteps*startDay:daySteps*endDay], 'r:', label = 'GW')
	plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['psi_w'][daySteps*startDay:daySteps*endDay], 'b-', label = 'Storage')
	plt.legend(['Soil','Base','Node', 'Leaf', 'GW', 'Storage'])
	plt.minorticks_on()
	plt.grid(visible=True, which='both', axis='both')
	
	# Save all figures to the same folder as resultsFile
	results_dir = os.path.dirname(resultsFile)
	print(results_dir)
	os.makedirs(results_dir, exist_ok=True)
	
	plot_suffix = f'_cw{plant_s[x]}_cs{s_s[x]}_cgw{gw_s[x]}'
	
	anp.savefig(os.path.join(results_dir, f'SoilMoisture{plot_suffix}.png'), dpi=300, bbox_inches='tight')
	fluxes.savefig(os.path.join(results_dir, f'Fluxes{plot_suffix}.png'), dpi=300, bbox_inches='tight')
	potents.savefig(os.path.join(results_dir, f'Potentials{plot_suffix}.png'), dpi=300, bbox_inches='tight')



