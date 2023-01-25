import numpy as np
import mesa_reader as mr
from model_io import load_orbital_state
import params

import sys

def main():
	# Check whether cur_time is outside the times spanned by the current MESA profiles

	cur_path=sys.argv[1]

	# load cur_time
	oh_finame="orbital_history.data"
	cur_time,_,_,_=load_orbital_state(oh_finame)

	# load times of MESA profiles
	if params.use_stored_profiles:
		# assume only local profiles exist
		sh_finame="{}/LOGS/history_full.data".format(cur_path)
	else:
		# assume we have all profiles
		# we should always be in-bounds until the simulation finishes (assuming t0 set properly)
		sh_finame="{}/LOGS/history.data".format(cur_path)
	h=mr.MesaData(sh_finame)

	# compare times
	if (h.star_age[0]<=cur_time) and (h.star_age[-1]>=cur_time):
		# in bounds
		sys.exit(0)
	else:
		# out of bounds
		sys.exit(1)

if __name__=="__main__":
	main()
