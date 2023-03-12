import sys

import numpy as np
import mesa_reader as mr

import params
from model_io import update_orbital_parameters, load_orbital_state
from unit_conversion import a_to_OmegaOrb

def main():
	# load command line arguments
	pind=sys.argv[1]
	cur_param_ind=int(sys.argv[2])
	cur_M=float(sys.argv[3])

	if pind==1:
		# Initialize orbital configuration in the GYRE inlist with the initial parameters from the params file
		update_orbital_parameters(params.OmegaOrb0[cur_param_ind], params.OmegaRot0[cur_param_ind], params.e0[cur_param_ind], params.gyre_inlist)
	else:
		# load the orbital state from the orbital history file
		oh_finame="orbital_history.data"
		cur_time,cur_a,cur_e,cur_OmegaRot=load_orbital_state(oh_finame)
		# convert semi-major axis to orbital frequency
		cur_OmegaOrb=a_to_OmegaOrb(cur_a, cur_M) # [cyc/day]

		# Update orbital parameters in the gyre_tides inlist
		update_orbital_parameters(cur_OmegaOrb, cur_OmegaRot, cur_e, params.gyre_inlist)

if __name__=="__main__":
	main()