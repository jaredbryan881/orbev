import numpy as np
import mesa_reader as mr
from model_io import load_orbital_state

import sys

def main():
	# Check whether cur_time is outside the times spanned by the current MESA profiles

	# load cur_time
	oh_finame="orbital_history.data"
	cur_time,_,_,_=load_orbital_state(oh_finame)

	# load times of MESA profiles
	sh_finame="./LOGS/history.data"
	h=mr.MesaData(sh_finame)

	# compare times
	if (h.star_age[0]<cur_time) and (h.star_age[-1]>cur_time):
		# out of bounds
		sys.exit(1)
	else:
		# not out of bounds
		sys.exit(0)

if __name__=="__main__":
	main()