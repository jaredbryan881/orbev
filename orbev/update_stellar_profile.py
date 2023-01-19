import numpy as np
import sys

from model_io import load_profile, save_profile, load_orbital_state
from interpolate_profile import select_profiles, get_interpolation_axis, lin_interp_2d, interpolate_single_quantity
import mesa_reader as mr

import params

def main():
	# Read stellar history file
	sh_finame="./LOGS/history_full.data"
	sh=mr.MesaData(sh_finame)

	# Read orbital configuration
	oh_finame="orbital_history.data"
	cur_time,cur_a,cur_e,cur_OmegaRot=load_orbital_state(oh_finame)

	# First, check whether we are close enough to a grid point to just load an exact stellar model
	# distance between two nearest profiles
	pnum1,pnum2=select_profiles(cur_time, sh.star_age, params.allowable_profiles, 2)
	pct=(cur_time-sh.star_age[pnum1])/(sh.star_age[pnum2]-sh.star_age[pnum1])
	assert (pct>=0 and pct<=1), "{} is not a valid range for interpolation. Must be between 0 and 1.".format(pct)
	print("Currently {}% between profile {} and {}".format(pct*100, pnum1+1, pnum2+1))
	# load profiles and their headers
	if pct<0.001:
		print("Loading profile {}".format(pnum1+1))
		p,header=load_profile("./LOGS/profile{}.data.GYRE".format(pnum1+1))
		header=np.array([[int(header[0]), header[1], header[2], header[3], int(header[4])]])
		save_profile(p,header)
	elif (1-pct)<0.001:
		print("Loading profile {}".format(pnum2+1))
		p,header=load_profile("./LOGS/profile{}.data.GYRE".format(pnum2+1))
		header=np.array([[int(header[0]), header[1], header[2], header[3], int(header[4])]])
		save_profile(p,header)
	else:
		# If we weren't close enough to a grid point, then we need to interpolate between the stellar models
		print("Interpolating {}% between profile {} and {}".format(pct*100, pnum1+1, pnum2+1))
		p1,header1=load_profile("./LOGS/profile{}.data.GYRE".format(pnum1+1))
		p2,header2=load_profile("./LOGS/profile{}.data.GYRE".format(pnum2+1))

		r1_interp, r2_interp=get_interpolation_axis(p1["r"], p2["r"])

		p_mid={}
		for key in p1.keys():
			if key=="ind":
				p_mid["ind"]=np.arange(1,len(r1_interp)+1)
			elif key=="r":
				p_mid["r"]=lin_interp_2d(r1_interp, r2_interp, pct)
			else:
				p_mid[key]=interpolate_single_quantity(r1_interp, p1, r2_interp, p2, pct, key)

		# interpolate header values
		M_mid=lin_interp_2d(header1[1], header2[1], pct)
		R_mid=lin_interp_2d(header1[2], header2[2], pct)
		L_mid=lin_interp_2d(header1[3], header2[3], pct)

		# write interpolated header and profile to a profile.data.GYRE file
		header_mid=np.array([[len(r1_interp), M_mid, R_mid, L_mid, int(101)]])
		save_profile(p_mid, header_mid)

	# interpolate and save current stellar moment of inertia
	Is=np.loadtxt("./LOGS/stellar_MOIs.txt")
	cur_I = lin_interp_2d(Is[pnum1], Is[pnum2], pct)
	np.savetxt("current_stellar_MOI.txt", np.array([cur_I]))

if __name__=="__main__":
	main()
