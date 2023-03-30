import numpy as np

###########
# Toggles #
###########
live_orbit=True # evolve the orbital configuration?
time_reversed=False # evolve time backward?
use_stored_profiles=False # assume we have all the profiles we need? or else generate them as we go.

#################################
# Initial Orbital Configuration #
#################################
data=np.loadtxt("params_ps_orot.txt")
e0=data[0]
OmegaOrb0=data[1]
OmegaRot0=data[2]
t0=None#data[3]*1e9
Ms=data[4]
Zs=data[5]

####################################
# Parameters for adaptive timestep #
####################################
max_dt = 1e6 # [yr]
base_dt = 1e4 # [yr]
e_eps = 1e-6 # max allowable fractional error in e
a_eps = 1e-6 # mac allowable fractional error in a
OmegaRot_eps = 1e-6 # max allowable fractional error in OmegaRot
safety_factor = 0.9 # Constant slightly <1 used to set adaptive timestep

#######################
# Paths and Filenames #
#######################
gyre_inlist="./gyre_tides.in"

# subset of profiles which can be selected from
allowable_profiles=np.arange(1,1006,1)

def main():
	import sys
	ip_ind=int(sys.argv[1])
	print("####################################")
	print("#### Initial Orbital Parameters ####")
	print("####################################")
	if t0 is None:
		print("t_0: ZAMS")
	else:
		print("t_0: {} Gyr".format(t0[ip_ind]/1e9))
	print("e_0: {}".format(e0[ip_ind]))
	print("Omega_orb0: {} cyc/day".format(OmegaOrb0[ip_ind]))
	print("Omega_rot0: {} cyc/day".format(OmegaRot0[ip_ind]))

if __name__=="__main__":
	main()
