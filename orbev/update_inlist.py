import sys

import numpy as np
import mesa_reader as mr

import params
from model_io import update_orbital_parameters, load_orbital_state, update_history
from unit_conversion import a_to_OmegaOrb

def main():
	# load command line arguments
	pind=int(sys.argv[1])
	rk_ind=int(sys.argv[2])
	cur_param_ind=int(sys.argv[3])
	cur_M=float(sys.argv[4])

	if rk_ind==0:
		# load the orbital state from the orbital history file
		cur_time,cur_a,cur_e,cur_OmegaRot,cur_dt=load_orbital_state("orbital_history.data")

		# TODO: move this update_history call somewhere else. Would be better if this script only updated the inlist
		# initialize another orbital history file just for the RKF substeps
		update_history(cur_time, cur_a, cur_e, cur_OmegaRot, cur_dt, foname="RKF_buffer.data")
	else:
		# calculate updated orbital parameters based on RKF step
		cur_time,cur_a,cur_e,cur_OmegaRot,cur_dt=load_orbital_state("RKF_buffer.data")
	
	# convert semi-major axis to orbital frequency
	cur_OmegaOrb=a_to_OmegaOrb(cur_a, cur_M) # [cyc/day]
	
	# Update orbital parameters in the gyre_tides inlist
	update_orbital_parameters(cur_OmegaOrb, cur_OmegaRot, cur_e, params.gyre_inlist)
		
if __name__=="__main__":
	main()