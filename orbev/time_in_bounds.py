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
		# assume we have every profile stored locally
		sh_finame="{}/star_ages.txt".format(foname)
		_,star_age=np.loadtxt(sh_finame)
	else:
		# take only times from current segment of profiles
		sh_finame="{}/history.data".format(foname)
		sh=mr.MesaData(sh_finame)
		star_age=sh.star_age

	if params.time_reversed:
		# compare times
		if (np.min(star_age)<=cur_time-cur_dt) and (np.max(star_age)>=(cur_time)):
			# in bounds
			sys.exit(0)
		else:
			# out of bounds
			print("Time ({} Myr) out of bounds ({}-{} Myr)".format(cur_time/1e6, star_age[0]/1e6, star_age[-1]/1e6))
			sys.exit(1)
	else:
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
