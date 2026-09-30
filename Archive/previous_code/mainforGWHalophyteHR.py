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
#from GWHalophyteHRmass_in import *
from GWHalophyteHR import *
import openpyxl
import time

duration = 30
weatherFile = "sample_data/Almeria2017JulyAugust.xlsx"
resultsLocation = "sample_output/SolverTrials"
resultsFile = resultsLocation+ '/November20.csv' # default value for the location where results are saved

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
sinit = 0.12 #Initial soil moisture
vwi = 0.86 #0.81 #Initial plant storage volume
gw_z_init = 3 #Initial depth to the groundwater table (m)
cs_init = 200 #200 #initial salt concentration in soil, (Mol/m^3)
cgw_init = 300 #(from TDS conversion of Armas data) ##initial salt concentration in groundwater (Mol/m^3)
cw_init = 270 #270 #initial salt concentration in plant storage (Mol/m^3)
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
print(results['hr_cum'][-1])
	
data = pd.DataFrame.from_dict(results)
datacsv = pd.DataFrame.from_dict(results)
datacsv.to_csv(resultsFile)
#data.to_pickle(resultsFile)


#Plot results
startDay = 0
endDay = duration
dispDuration = endDay-startDay
daySteps = 60//timestepM*24
timevec = np.linspace(0,duration,duration*daySteps)
timevecHr = np.linspace(0,duration*24,duration*daySteps)


anp = plt.figure()
gph_title = 'Soil moisture: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init)
plt.title(gph_title)
plt.xlabel("time (days)")
plt.ylabel("s (-)")
plt.plot(timevec, results['s'][daySteps*startDay:daySteps*endDay])
for i in range(0, endDay):
 	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.xlim(0, dispDuration)
plt.ylim(0, 0.3)
plt.xticks([0.,6.,12.,18.,24.,30.])
plt.grid(true, axis='y')
plt.legend()
#anp.savefig(resultsLocation +'/soil moisture.png')
anp.show()

fluxes = plt.figure()
fluxes.suptitle('Fluxes: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
plt.xlabel("t (days)")
plt.ylabel("Flux (um/s)")
plt.ylim([min(results['qs'])-0.01,max(results['ev'])+0.01])
plt.xlim([0,endDay])
print(timevec)
plt.plot(timevec, results['qs'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
#plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], results['qs'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
plt.plot(timevec, results['qgw'][daySteps*startDay:daySteps*endDay], 'k-', label = 'GW')
plt.plot(timevec, results['ev'][daySteps*startDay:daySteps*endDay], 'b:', label = 'ET')
plt.plot(timevec, results['qw'][daySteps*startDay:daySteps*endDay], 'g:', label = 'Storage')
plt.legend(['Soil','GW','EV', 'Storage'], loc='upper right')
plt.minorticks_on()
plt.grid(b ='true', which = 'both', axis = 'both')
#fluxes.savefig(resultsLocation + '/fluxes.png')

potents = plt.figure()
potents.suptitle('Potentials: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
plt.xlabel('time (days)')
plt.ylabel('Potential (MPa)')
plt.xlim([0, endDay])
plt.ylim([-4.5,0.5])
plt.plot(timevec, results['psi_s'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
plt.plot(timevec, results['psi_b'][daySteps*startDay:daySteps*endDay], 'k-', label = 'Base')
plt.plot(timevec, results['psi_x'][daySteps*startDay:daySteps*endDay], 'b:', label = 'Node')
plt.plot(timevec, results['psi_l'][daySteps*startDay:daySteps*endDay], 'g-', label = 'Leaf')
plt.plot(timevec, results['psi_gw'][daySteps*startDay:daySteps*endDay], 'r:', label = 'GW')
plt.plot(timevec, results['psi_w'][daySteps*startDay:daySteps*endDay], 'c:', label = 'PWS')
plt.legend(['Soil','Base','Node', 'Leaf', 'GW', 'PWS'])
plt.minorticks_on()
plt.grid(b ='true', which = 'both', axis = 'both')
#potents.savefig(resultsLocation +'/potentials.png')

##### Potential Plots of Plant, Soil, Subsoil
total_potent = plt.figure()
#total_potent.suptitle('Total Potential: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
plt.xlabel('time (days)')
plt.ylabel('Potential (MPa)')
plt.xlim([0, endDay])
plt.ylim([-4.5,0])
plt.plot(timevec, results['psi_s'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
plt.plot(timevec, results['psi_gw'][daySteps*startDay:daySteps*endDay], 'r-', label = 'GW')
plt.plot(timevec, results['psi_w'][daySteps*startDay:daySteps*endDay], 'k-', label = 'Storage')
plt.legend(['Soil','GW', 'Storage'])
plt.minorticks_on()
plt.grid(b ='true', which = 'both',axis = 'both')
#total_potent.savefig(resultsLocation + '/total potentials.png')

osm_potent = plt.figure()
osm_potent.suptitle('Osmotic Potential: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
plt.xlabel('time (days)')
plt.ylabel('Potential (MPa)')
plt.xlim([0, endDay])
plt.ylim([-3,0.5])
plt.plot(timevec, results['psi_s_osm'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
plt.plot(timevec, results['psi_gw_osm'][daySteps*startDay:daySteps*endDay], 'r:', label = 'GW')
plt.plot(timevec, results['psi_w_osm'][daySteps*startDay:daySteps*endDay], 'c:', label = 'PWS')
plt.legend(['Soil','GW', 'PWS'])
plt.minorticks_on()
plt.grid(b ='true', which = 'both',axis = 'both')
#osm_potent.savefig(resultsLocation + '/osmotic potentials.png')

mat_potent = plt.figure()
mat_potent.suptitle('Matric/Turgor Potential: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
plt.xlabel('time (days)')
plt.ylabel('Potential (MPa)')
plt.xlim([0, endDay])
plt.ylim([-4.5, 2])
plt.plot(timevec, results['psi_s_mat'][daySteps*startDay:daySteps*endDay], 'm-', label = 'Soil')
plt.plot(timevec, results['psi_gw_mat'][daySteps*startDay:daySteps*endDay], 'r:', label = 'GW')
plt.plot(timevec, results['psi_w_turgor'][daySteps*startDay:daySteps*endDay], 'c:', label = 'PWS')
plt.legend(['Soil','GW', 'PWS'])
plt.minorticks_on()
plt.grid(b ='true', which = 'both',axis = 'both')
#mat_potent.savefig(resultsLocation +'/matric potentials.png')

hr_cum = plt.figure()
hr_cum.suptitle('Cumulative HR: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
plt.xlabel('time (days)')
plt.ylabel('Cumulative HR (mm)')
plt.xlim([0, endDay])
#plt.ylim([-3,0.5])
plt.plot(timevec, results['hr_cum'][daySteps*startDay:daySteps*endDay], 'm-')
plt.minorticks_on()
plt.grid(b ='true', which = 'both',axis = 'both')
#hr_cum.savefig(resultsLocation + '/hr cumulative.png')

hr_matric = plt.figure()
hr_matric.suptitle('Cumulative HR vs Soil Matric Potential: cw {} cs {} cgw {}'.format(cw_init, cs_init, cgw_init))
plt.xlabel('Matric Potential (-MPa)')
plt.ylabel('Cumulative HR (mm)')
#plt.xlim([0, endDay])
#plt.ylim([-3,0.5])
plt.plot([i*(-1) for i in results['psi_s_mat'][daySteps*startDay:daySteps*endDay]], results['hr_cum'][daySteps*startDay:daySteps*endDay], 'k-')
plt.minorticks_on()
plt.grid(b ='true', which = 'both',axis = 'both')



# =============================================================================
# """3d plots"""
# phi_i = 148 #Value from Starry2014 Weather Set
# ta_i = 21 #Value from Starry2014 Weather Set
# ta_i_K = ta_i+273 #Atmospheric temperature in K
# rh_i = 53 #Value from Starry2014 Weather Set
# psat_i = A_SAT*np.exp((B_SAT*(ta_i))/(C_SAT + ta_i))
# qaInp_i = 0.622*rh_i/100.*psat_i/P_ATM
# cm_i = 350
# vw_i = 0.030511169 #Value at 15th time step
# psi_w_i = -0.212625941
# cw_i = 150.1206692
# psi_l_env = np.linspace(-3, 0, 200)
# tl_env = np.linspace(287, 295, 200)
# Psi_l_env, Tl_env = np.meshgrid(psi_l_env, tl_env)
# 
# 
# #Ev calculated for envelope of psi_l, tl values
# Ev_env = []
# for i in range(len(psi_l_env)):
# 	start = time.time()
# 	Ev_env.append([])
# 	for j in range(len(tl_env)):
# 		Ev_env[i].append(hydro.evf(photo, phi_i, ta_i_K, psi_l_env[i], qaInp_i, tl_env[j], photo.cm, hydro.lai, 1))
# 	end = time.time()
# 	#print(end-start)
# 
# water_flux = []
# for i in range(len(psi_l_env)):
# 	start = time.time()
# 	water_flux.append([])
# 	for j in range(len(tl_env)):
# 		water_flux[i].append(hydro.qbx(hydro.gp, hydro.psi_x(hydro.evf(photo, phi_i, ta_i_K, psi_l_env[i], qaInp_i, tl_env[j], photo.cm, hydro.lai, 1),psi_l_env[i],hydro.gp,hydro.lai), hydro.psi_b(vw_i,hydro.evf(photo, phi_i, ta_i_K, psi_l_env[i], qaInp_i, tl_env[j], photo.cm, hydro.lai, 1), psi_l_env[i],psi_w_i,hydro.gp,hydro.gwf(psi_w_i),hydro.lai,cw_i,soil),hydro.lai) + hydro.qwf(vw_i,hydro.evf(photo, phi_i, ta_i_K, psi_l_env[i], qaInp_i, tl_env[j], photo.cm, hydro.lai, 1),hydro.gp,psi_l_env[i],hydro.lai, cw_i,dt))
# 
# 
# 
# ev3d = plt.figure()
# ax3d_1 = plt.axes(projection='3d')
# #ax3d_1.contour3D(psi_l_env, tl_env, np.array(Ev_env), 50, cmap='binary')
# #ax3d_1.plot_surface(Psi_l_env, Tl_env, np.array(Ev_env), rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax3d_1.plot_wireframe(Psi_l_env, Tl_env, np.array(Ev_env), color='black')
# ax3d_1.plot_wireframe(Psi_l_env, Tl_env, np.array(water_flux), color='red')
# ax3d_1.set_xlabel('psi_l (MPa)')
# ax3d_1.set_ylabel('tl (K)')
# ax3d_1.set_zlabel('ev/hydro_flux (um/s)')
# ax3d_1.set_zlim3d(0,0.4)
# ax3d_1.view_init(45,160)
# =============================================================================


#print(hydro.evf(photo, phi_i, ta_i_K, -0.847017435501732, qaInp_i, 293.80209, photo.cm, hydro.lai, 1)) Sanity check




