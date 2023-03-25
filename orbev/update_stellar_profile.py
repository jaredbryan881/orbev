import numpy as np
import sys

from model_io import load_profile, save_profile, load_orbital_state
from interpolate_profile import select_profiles, get_interpolation_axis, lin_interp_2d, interpolate_single_quantity
import mesa_reader as mr

import pickle as pkl

import params

from RK45 import RK45_tableau

def main():
	base_profile_dir=sys.argv[1]
	rk_ind=int(sys.argv[2])

	# Read stellar history file
	# don't confuse this with history_full.data, whose indices would not relate to pnums
	sh_finame="{}/history.data".format(base_profile_dir)
	sh=mr.MesaData(sh_finame)

	# Read orbital configuration
	oh_finame="orbital_history.data"
	cur_time,_,_,_,cur_dt=load_orbital_state(oh_finame)

	# update the current time based on where we are in the RKF step
	tab=RK45_tableau()
	cur_time = cur_time + tab.c[rk_ind]*cur_dt

	# First, check whether we are close enough to a grid point to just load an exact stellar model
	# distance between two nearest profiles
	pnum1,pnum2=select_profiles(cur_time, sh.star_age, params.allowable_profiles, 2)
	pct=(cur_time-sh.star_age[pnum1])/(sh.star_age[pnum2]-sh.star_age[pnum1])
	assert (pct>=0 and pct<=1), "{} is not a valid range for interpolation. Must be between 0 and 1.".format(pct)

	# load profiles and their headers
	with open("{}/profiles.pkl".format(base_profile_dir), "rb") as f:
		headers,profiles=pkl.load(f)

	# select the current profiles
	if pct<0.001:
		print("Loading profile {}".format(pnum1+1))
		profile=profiles[pnum1]
		header=headers[pnum1]
		header=np.array([int(header[0]), header[1], header[2], header[3], int(header[4])])
		save_profile(profile,header)
	elif (1-pct)<0.001:
		print("Loading profile {}".format(pnum2+1))
		profile=profiles[pnum2]
		header=headers[pnum2]
		header=np.array([int(header[0]), header[1], header[2], header[3], int(header[4])])
		save_profile(profile,header)
	else:
		# If we weren't close enough to a grid point, then we need to interpolate between the stellar models
		print("Interpolating {}% between profile {} and {}".format(pct*100, pnum1+1, pnum2+1))
		profile1=profiles[pnum1]
		profile2=profiles[pnum2]

		header1=headers[pnum1]
		header1=np.array([int(header1[0]), header1[1], header1[2], header1[3], int(header1[4])])
		header2=headers[pnum2]
		header2=np.array([int(header2[0]), header2[1], header2[2], header2[3], int(header2[4])])

		r1_interp, r2_interp=get_interpolation_axis(profile1["r"], profile2["r"])

		p_mid={}
		for key in profile1.keys():
			if key=="ind":
				p_mid["ind"]=np.arange(1,len(r1_interp)+1)
			elif key=="r":
				p_mid["r"]=lin_interp_2d(r1_interp, r2_interp, pct)
			else:
				p_mid[key]=interpolate_single_quantity(r1_interp, profile1, r2_interp, profile2, pct, key)

		# interpolate header values
		M_mid=lin_interp_2d(header1[1], header2[1], pct)
		R_mid=lin_interp_2d(header1[2], header2[2], pct)
		L_mid=lin_interp_2d(header1[3], header2[3], pct)

		# write interpolated header and profile to a profile.data.GYRE file
		header_mid=np.array([[len(r1_interp), M_mid, R_mid, L_mid, int(101)]])
		save_profile(p_mid, header_mid)

	# interpolate and save current stellar moment of inertia
	Is=np.loadtxt("{}/stellar_MOIs.txt".format(base_profile_dir))
	cur_I = lin_interp_2d(Is[pnum1], Is[pnum2], pct)
	np.savetxt("current_stellar_MOI.txt", np.array([cur_I]))

if __name__=="__main__":
	main()
