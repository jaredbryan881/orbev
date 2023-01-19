# Take the current orbital configuration and maximum allowable time/parameter steps
# Return the timestep and updated orbital parameters

import sys
import os

import numpy as np
import h5py
import pandas as pd

import mesa_reader as mr

import params
from model_io import load_profile, load_stellar_state, load_orbital_state, update_orbital_parameters, update_history
from unit_conversion import freq_scale, a_to_OmegaOrb, OmegaOrb_to_a
from calculate_Is import MOI

# constants
G = 6.67430e-11  # [m3 kg-1 s-2]
Rsun = 695700000 # [m]
Msun = 1.989e30  # [kg]
au = 1.496e+11   # [m]

def main():
	# load command line arguments
	pind=int(sys.argv[1])
	cur_dir=sys.argv[2]
	cur_param_ind=int(sys.argv[3])

	# Read stellar history file
	sh_finame="./LOGS/history_full.data"
	sh=mr.MesaData(sh_finame)

	if pind==1:
		# Initialize orbital configuration in the GYRE inlist
		update_orbital_parameters(params.OmegaOrb0[cur_param_ind], params.OmegaRot0[cur_param_ind], params.e0[cur_param_ind], params.gyre_inlist)
		# Initialize orbital configuration history file
		# little hack to get around not having the orbital history file yet: 
		# just interpret the initial M from the cur_dir string as a mass in Msun
		# it's pretty ugly though
		cur_M=float(cur_dir.split('M')[1].split("_")[0])
		a0=OmegaOrb_to_a(params.OmegaOrb0[cur_param_ind], cur_M)
		# initialize history file
		update_history(params.t0, a0, params.e0[cur_param_ind], params.OmegaRot0[cur_param_ind])
		return

	# get stellar state from the current profile.data.GYRE
	cur_R,cur_M=load_stellar_state("profile_cur.data.GYRE") # [cm], [g]
	cur_R=cur_R/100/Rsun  # [Rsun]
	cur_M=cur_M/1000/Msun # [Msun]

	# Read orbital configuration
	oh_finame="orbital_history.data"
	cur_time,cur_a,cur_e,cur_OmegaRot=load_orbital_state(oh_finame)
	cur_OmegaOrb=a_to_OmegaOrb(cur_a, cur_M)

	# Define dimensionalizing constant
	fs = freq_scale(cur_M*Msun, cur_R*Rsun)

	# read tides output
	print("Loading tide_orbit")
	with h5py.File("./profile{}/tide_orbit.h5".format(pind-1), "r") as hf:
		# orbital parameters
		Omega_orb=hf["Omega_orb"][0]*fs

		a=OmegaOrb_to_a(Omega_orb, cur_M)

		e=hf["e"][0]

		# orbital evolution rates
		edot=hf["e_dot"][0]*fs*365 # 1/yr

		adot=hf["a_dot"][0]*fs*365*cur_R*(Rsun/au) # au/yr

		Jdot=hf["J_dot"][0]*fs*365 # 1/yr
		Jdot*=(cur_R**(1/2) * ((G*(Rsun**3)/Msun)**(1/2))*86400 * cur_M**(3/2)) # Msun*Rsun**2*cyc/day/yr
		profile, header = load_profile("./profile_cur.data.GYRE")

		if os.path.exists("current_stellar_MOI.txt"):
			I = np.loadtxt("current_stellar_MOI.txt")
		else:
			I = MOI(profile["M"], profile["r"]) # Msun*Rsun^2
		OmegaRotdot = Jdot/I # cyc/day/yr

	# calculate the timestep and update orbital parameters
	if params.live_orbit:
		# calculate the timestep
		dt = np.min([np.abs(params.max_de/edot), np.abs(params.max_da/adot), np.abs(params.max_dOmegaRot/OmegaRotdot), params.max_dt])
		# step backward in time if time-reversed simulation
		if params.time_reversed:
			dt=-dt
		new_time=cur_time+dt

		# update orbital parameters
		new_e = e + dt*edot
		new_a = a + dt*adot
		new_OmegaOrb = a_to_OmegaOrb(new_a, cur_M)
		new_OmegaRot = cur_OmegaRot + dt*OmegaRotdot
	else:
		# take fixed timestep
		dt = params.max_dt
		# step backward in time if time-reversed simulation
		if params.time_reversed:
			dt=-dt
		new_time = cur_time+dt

		# update orbital parameters
		new_e = params.e0[cur_param_ind]
		new_a = OmegaOrb_to_a(params.OmegaOrb0[cur_param_ind], cur_M)
		new_OmegaOrb = params.OmegaOrb0[cur_param_ind]
		new_OmegaRot = params.OmegaRot0[cur_param_ind]

	print("Time: {} Myr".format(cur_time/1e6))
	print("dt = {} Myr".format(dt/1e6))
	print("edot = {} 1/yr; de={}".format(edot, dt*edot))
	print("adot = {} au/yr; da={}".format(adot, dt*adot))
	print("OmegaRotdot = {} cyc/day/yr; dOmegaRot={}".format(OmegaRotdot, dt*OmegaRotdot))

	# update orbital_history.data
	update_history(new_time, new_a, new_e, new_OmegaRot)
	# update parameters in gyre_orbit.in
	update_orbital_parameters(new_OmegaOrb, new_OmegaRot, new_e, params.gyre_inlist)

if __name__=="__main__":
	main()