import sys

import numpy as np
import mesa_reader as mr

import params
from model_io import update_orbital_parameters, update_history
from unit_conversion import OmegaOrb_to_a

def main():
	# load command line arguments
	cur_path=sys.argv[1]
	cur_dir=sys.argv[2]
	cur_param_ind=int(sys.argv[3])

	# Initialize orbital configuration history file
	# Read stellar history file
	sh_finame="{}/LOGS/history_full.data".format(cur_path)
	sh=mr.MesaData(sh_finame)
	# little hack to get around not having the orbital history file yet: 
	# just interpret the initial M from the cur_dir string as a mass in Msun
	# it's pretty ugly though
	cur_M=float(cur_dir.split('M')[1].split("_")[0])
	a0=OmegaOrb_to_a(params.OmegaOrb0[cur_param_ind], cur_M)
	# initialize history file
	if params.t0 is None:
		t0=np.ceil(sh.star_age[0])
	else:
		t0=params.t0[cur_param_ind]
	update_history(t0, a0, params.e0[cur_param_ind], params.OmegaRot0[cur_param_ind], params.base_dt)

if __name__=="__main__":
	main()