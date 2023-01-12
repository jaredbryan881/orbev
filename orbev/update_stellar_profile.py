import numpy as np
import sys

from model_io import load_profile, save_profile, load_orbital_state
from interpolate_profile import select_profiles, resample_radial_axes, resample_profiles, interp_2d
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
	
	# If we weren't close enough to a grid point, then we need to interpolate between the stellar models
	# find the N nearest stellar models in time
	selected_pnums = select_profiles(cur_time, sh.star_age, params.allowable_profiles, params.n_profiles)

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
		print("Currently {}% between profile {} and {}. Performing cubic spline interpolation with {} profiles.".format(pct*100, pnum1, pnum2, params.n_profiles))

		# load the selected profiles and their headers
		profiles=[]
		headers=[]
		for pnum in selected_pnums:
			p, h = load_profile(base_sh_finame+"profile{}.data.GYRE".format(pnum+1))
			profiles.append(p)
			headers.append(h)

		# get each profiles in selected_profiles to have the same number of samples
		# done by finding the profile with the most samples and cascading Monge maps to couple each profile
		rs_interp = resample_radial_axes([p["r"] for p in profiles])
		profiles_interp = resample_profiles(rs_interp, profiles)

		p_mid={}
		for key in profiles_interp[0].keys():
			if key=="ind":
				p_mid["ind"]=np.arange(1,len(rs_interp[0])+1)
			elif key=="r":
				p_mid["r"]=interp_2d(rs_interp, cur_time, sh.star_age[selected_pnums], params.spline_order)
			else:
				p_mid[key]=interp_2d([p[key] for p in profiles_interp], cur_time, sh.star_age[selected_pnums], params.spline_order)

		# interpolate header values
		M_mid=interp_2d([h[1] for h in headers], cur_time, sh.star_age[selected_pnums], params.spline_order)
		R_mid=interp_2d([h[2] for h in headers], cur_time, sh.star_age[selected_pnums], params.spline_order)
		L_mid=interp_2d([h[3] for h in headers], cur_time, sh.star_age[selected_pnums], params.spline_order)

		# write interpolated header and profile to a profile.data.GYRE file
		header_mid=np.array([[len(rs_interp[0]), M_mid, R_mid, L_mid, int(101)]])
		save_profile(p_mid, header_mid)

	# interpolate and save current stellar moment of inertia
	Is=np.loadtxt(base_sh_finame+"stellar_MOIs.txt")
	cur_I = interp_2d(Is[selected_pnums], cur_time, sh.star_age[selected_pnums], params.spline_order)
	np.savetxt("current_stellar_MOI.txt", np.array([cur_I]))

if __name__=="__main__":
	main()
