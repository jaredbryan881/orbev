import numpy as np
import mesa_reader as mr
from model_io import load_orbital_state
import params

import sys

def main():
	# Check whether cur_time is outside the times spanned by the current MESA profiles

	foname=sys.argv[1]

	# load cur_time
	oh_finame="orbital_history.data"
	cur_time,_,_,_,cur_dt=load_orbital_state(oh_finame)

	# load times of MESA profiles
	if params.use_stored_profiles:
		# assume only local profiles exist
		sh_finame="{}/star_ages.txt".format(foname)
		star_age=np.loadtxt(sh_finame)
	else:
		# assume we have all profiles
		# we should always be in-bounds until the simulation finishes (assuming t0 set properly)
		sh_finame="{}/history.data".format(foname)
		h=mr.MesaData(sh_finame)
		star_age=h.star_age

	# compare times
	if (np.min(star_age)<=cur_time) and (np.max(star_age)>=(cur_time+cur_dt)):
		# in bounds
		sys.exit(0)
	else:
		# out of bounds
		print("Time ({} Myr) out of bounds ({}-{} Myr)".format(cur_time/1e6, star_age[0]/1e6, star_age[-1]/1e6))
		sys.exit(1)

if __name__=="__main__":
	main()
