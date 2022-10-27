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

	# Read stellar history file
	sh_finame="/home/jared/MIT/astero/mesa_HATP2/live_planet/{}/LOGS/history.data".format(cur_dir)
	sh=mr.MesaData(sh_finame)

	# get stellar state from the current profile.data.GYRE
	cur_R,cur_M=load_stellar_state("profile_cur.data.GYRE") # [cm], [g]
	cur_R=cur_R/100/Rsun  # [Rsun]
	cur_M=cur_M/1000/Msun # [Msun]

	if pind==1:
		# Initialize orbital configuration in the GYRE inlist
		update_orbital_parameters(params.OmegaOrb0, params.OmegaRot0, params.e0, params.finame)
		# Initialize orbital configuration history file
		a0=OmegaOrb_to_a(params.OmegaOrb0, cur_M)
		update_history(sh.star_age[0], a0, params.e0, params.OmegaRot0)
		return

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
		Omega_orb=hf["Omega_orb"][0]
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

	# calculate the timestep
	dt = np.min([np.abs(params.max_de/edot), np.abs(params.max_da/adot), np.abs(params.max_dOmegaRot/OmegaRotdot), params.max_dt])
	new_time=cur_time+dt

	print("Time: {} Myr".format(cur_time/1e6))
	print("dt = {} Myr".format(dt/1e6))
	print("edot = {} 1/yr".format(edot))
	print("adot = {} au/yr".format(adot))
	print("OmegaRotdot = {} cyc/day/yr".format(OmegaRotdot))

	# update e
	new_e = e + dt*edot
	# update a (and thus OmegaOrb)
	new_a = a + dt*adot
	new_OmegaOrb = a_to_OmegaOrb(new_a, cur_M)
	# update OmegaRot
	new_OmegaRot = cur_OmegaRot + dt*OmegaRotdot

	# update orbital_history.data
	update_history(new_time, new_a, new_e, new_OmegaRot)
	# update parameters in gyre_orbit.in
	update_orbital_parameters(new_OmegaOrb, new_OmegaRot, new_e, params.finame)

if __name__=="__main__":
	main()