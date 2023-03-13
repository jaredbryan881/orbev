# Take the current orbital configuration and maximum allowable time/parameter steps
# Return the timestep and updated orbital parameters

import sys
import os

import numpy as np
import h5py
import pandas as pd

import mesa_reader as mr

import params
from model_io import load_profile, load_stellar_state, load_orbital_state, update_history
from unit_conversion import freq_scale, a_to_OmegaOrb, OmegaOrb_to_a
from calculate_Is import MOI

from RK45 import RK45_tableau

# constants
G = 6.67430e-11  # [m3 kg-1 s-2]
Rsun = 695700000 # [m]
Msun = 1.989e30  # [kg]
au = 1.496e+11   # [m]

def main():
	# load command line arguments
	cur_path=sys.argv[1]
	cur_dir=sys.argv[2]
	cur_param_ind=int(sys.argv[3])
	orbev_fodir=sys.argv[4]
	rk_ind=int(sys.argv[5])

	# get stellar state from the current profile.data.GYRE
	cur_R,cur_M=load_stellar_state("profile_cur.data.GYRE") # [cm], [g]
	cur_R=cur_R/100/Rsun  # [Rsun]
	cur_M=cur_M/1000/Msun # [Msun]

	# Read orbital configuration
	cur_time,cur_a,cur_e,cur_OmegaRot,cur_dt=load_orbital_state("RKF_buffer.data")
	cur_OmegaOrb=a_to_OmegaOrb(cur_a, cur_M)

	# Define dimensionalizing constant
	fs = freq_scale(cur_M*Msun, cur_R*Rsun)

	# read tides output
	# [edot,adot,OmegaRotdot] make up the k vector for RK4(5)
	print("Loading tide_orbit")
	with h5py.File("{}/tidal_response_RKF_step.h5".format(orbev_fodir), "r") as hf:
		# orbital evolution rates
		edot=hf["e_dot"][:,0]*fs*365 # 1/yr

		adot=hf["a_dot"][:,0]*fs*365*cur_R*(Rsun/au) # au/yr

		Jdot=hf["J_dot"][:,0]*fs*365 # 1/yr
		Jdot*=(cur_R**(1/2) * ((G*(Rsun**3)/Msun)**(1/2))*86400 * cur_M**(3/2)) # Msun*Rsun**2*cyc/day/yr
		profile, header = load_profile("./profile_cur.data.GYRE")

		if os.path.exists("current_stellar_MOI.txt"):
			I = np.loadtxt("current_stellar_MOI.txt")
		else:
			I = MOI(profile["M"], profile["r"]) # Msun*Rsun^2
		OmegaRotdot = Jdot/I # cyc/day/yr

	# get an RK4(5) tableau to get the intermediate function evalution step size
	tab=RK45_tableau()

	# zero-pad the orbital evolution rates to make dot product possible
	edot=np.hstack([edot, np.zeros(5-rk_ind)])
	adot=np.hstack([adot, np.zeros(5-rk_ind)])
	OmegaRotdot=np.hstack([OmegaRotdot, np.zeros(5-rk_ind)])

	# calculate the timestep and update orbital parameters
	if params.live_orbit:
		# calculate the timestep
		# step backward in time if time-reversed simulation
		if params.time_reversed:
			cur_dt=-cur_dt
		new_time=cur_time+tab.c[rk_ind]*cur_dt

		# update orbital parameters
		new_e = cur_e + np.dot(tab.a[rk_ind,:],edot[:-1])*cur_dt
		new_a = cur_a + np.dot(tab.a[rk_ind,:],adot[:-1])*cur_dt
		new_OmegaOrb = a_to_OmegaOrb(new_a, cur_M)
		new_OmegaRot = cur_OmegaRot + np.dot(tab.a[rk_ind,:],OmegaRotdot[:-1])*cur_dt
	else:
		# TODO: update how we do fixed timesteps
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

	# update the RKF buffer
	update_history(new_time, new_a, new_e, new_OmegaRot, cur_dt, foname="RKF_buffer.data")

	# finally use the intermediate orbital evolution rates to make an update to orbital_history.data
	if rk_ind==5:
		retry_flag=False
		new_time = cur_time + cur_dt

		new_e = cur_e + cur_dt*np.dot(tab.b[0,:], edot)
		new_e_star = cur_e + cur_dt*np.dot(tab.b[1,:], edot)
		e_Delta1 = np.abs(new_e-new_e_star)
		e_Delta0 = params.e_eps*cur_e
		if e_Delta0>=e_Delta1:
			new_dt_e = params.safety_factor*cur_dt*(e_Delta0/e_Delta1)**(1/5)
		else:
			new_dt_e = params.safety_factor*cur_dt*(e_Delta0/e_Delta1)**(1/4)
			retry_flag=True

		new_a = cur_a + cur_dt*np.dot(tab.b[0,:], adot)
		new_a_star = cur_a + cur_dt*np.dot(tab.b[1,:], adot)
		a_Delta1 = np.abs(new_a-new_a_star)
		a_Delta0 = params.a_eps*cur_a
		if a_Delta0>=e_Delta1:
			new_dt_a = params.safety_factor*cur_dt*(a_Delta0/a_Delta1)**(1/5)
		else:
			new_dt_a = params.safety_factor*cur_dt*(a_Delta0/a_Delta1)**(1/4)
			retry_flag=True

		new_OmegaRot = cur_OmegaRot + cur_dt*np.dot(tab.b[0,:], OmegaRotdot)
		new_OmegaRot_star = cur_OmegaRot + cur_dt*np.dot(tab.b[1,:], OmegaRotdot)
		OmegaRot_Delta1 = np.abs(new_OmegaRot-new_OmegaRot_star)
		OmegaRot_Delta0 = params.OmegaRot_eps*cur_OmegaRot
		if OmegaRot_Delta0>=OmegaRot_Delta1:
			new_dt_OmegaRot = params.safety_factor*cur_dt*(OmegaRot_Delta0/OmegaRot_Delta1)**(1/5)
		else:
			new_dt_OmegaRot = params.safety_factor*cur_dt*(OmegaRot_Delta0/OmegaRot_Delta1)**(1/4)
			retry_flag=True

		new_dt = np.min([params.max_dt, new_dt_e, new_dt_a, new_dt_OmegaRot])
		print("New timestep is {}".format(new_dt))

		if retry_flag:
			print("Error too high, lowering timestep and trying again.")
			update_history(cur_time, cur_a, cur_e, cur_OmegaRot, new_dt, foname="orbital_history.data")
		else:
			update_history(new_time, new_a, new_e, new_OmegaRot, new_dt, foname="orbital_history.data")

if __name__=="__main__":
	main()