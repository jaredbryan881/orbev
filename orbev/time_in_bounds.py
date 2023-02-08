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
	cur_time,_,_,_=load_orbital_state(oh_finame)

	# load times of MESA profiles
	if params.use_stored_profiles:
		# assume only local profiles exist
		sh_finame="{}/history_full.data".format(foname)
	else:
		# assume we have all profiles
		# we should always be in-bounds until the simulation finishes (assuming t0 set properly)
		sh_finame="{}/history.data".format(foname)
	h=mr.MesaData(sh_finame)

	# compare times
	if (np.min(h.star_age)<=cur_time) and (np.max(h.star_age)>=cur_time):
		# in bounds
		sys.exit(0)
	else:
		# out of bounds
		print("Time ({} Myr) out of bounds ({}-{} Myr)".format(cur_time/1e6, h.star_age[0]/1e6, h.star_age[-1]/1e6))
		sys.exit(1)

if __name__=="__main__":
	main()
