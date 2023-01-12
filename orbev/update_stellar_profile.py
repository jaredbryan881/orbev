import numpy as np
import sys

from model_io import load_profile, save_profile, load_orbital_state
from interpolate_profile import select_profiles, get_interpolation_axis, lin_interp_2d, interpolate_single_quantity
import mesa_reader as mr

import params

def main():
	# load command line arguments
	pind=int(sys.argv[1])
	cur_dir=sys.argv[2]
	ip=int(sys.argv[3])

	# Read stellar history file
	base_sh_finame="{}/{}/LOGS/".format(params.mesa_diname, cur_dir)
	sh_finame=base_sh_finame+"history.data"
	sh=mr.MesaData(sh_finame)

	# Read orbital configuration
	oh_finame="orbital_history.data"
	cur_time,cur_a,cur_e,cur_OmegaRot=load_orbital_state(oh_finame)

	if pind==1:
		# just use the first profile copied in the shell script
		print("Starting with profile {}".format(ip))
		return

	# First, check whether we are close enough to a grid point to just load an exact stellar model
	# distance between two nearest profiles
	pnum1,pnum2=select_profiles(cur_time, sh.star_age, params.allowable_profiles, 2)
	pct=(cur_time-sh.star_age[pnum1])/(sh.star_age[pnum2]-sh.star_age[pnum1])
	
	# load profiles and their headers
	if pct<0.001:
		print("Loading profile {}".format(pnum1))
		p,header=load_profile(base_sh_finame+"profile{}.data.GYRE".format(pnum1))
		header=np.array([[int(header[0]), header[1], header[2], header[3], int(header[4])]])
		save_profile(p,header)
	elif (1-pct)<0.001:
		print("Loading profile {}".format(pnum2+1))
		p,header=load_profile(base_sh_finame+"profile{}.data.GYRE".format(pnum2))
		header=np.array([[int(header[0]), header[1], header[2], header[3], int(header[4])]])
		save_profile(p,header)
	else:
		# If we weren't close enough to a grid point, then we need to interpolate between the stellar models
		print("Interpolating {}% between profile {} and {}".format(pct*100, pnum1, pnum2))
		p1,header1=load_profile(base_sh_finame+"profile{}.data.GYRE".format(pnum1))
		p2,header2=load_profile(base_sh_finame+"profile{}.data.GYRE".format(pnum2))

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
	Is=np.loadtxt(base_sh_finame+"stellar_MOIs.txt")
	cur_I = lin_interp_2d(Is[pnum1], Is[pnum2], pct)
	np.savetxt("current_stellar_MOI.txt", np.array([cur_I]))

if __name__=="__main__":
	main()
