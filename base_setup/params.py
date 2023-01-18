import numpy as np

###########
# Toggles #
###########
live_orbit=False # evolve the orbital configuration?
time_reversed=False # evolve time backward?
store_all_profiles=False # store every single profile, or only those close to the current time?

#################################
# Initial Orbital Configuration #
#################################
t0=5e8 # [yr]
e0=np.array([0.5]*50) # []
OmegaOrb0=np.array([0.1774]*50) # [cyc/day]
OmegaRot0=np.array([0.1774]*50)*4.8 # [cyc/day]

#################
# Update limits #
#################
max_dt=1e4         # [yr]
max_de=1e-5        # []
max_da=1e-5        # [au]
max_dOmegaRot=1e-5 # [cyc/day]

#######################
# Paths and Filenames #
#######################
gyre_inlist="./gyre_tides.in"
mesa_diname="/home/jared/MIT/astero/mesa_HATP2/live_planet"

# subset of profiles which can be selected from
allowable_profiles=np.arange(1,1001,10)
n_profiles=50

def main():
	import sys
	ip_ind=int(sys.argv[1])
	print("####################################")
	print("#### Initial Orbital Parameters ####")
	print("####################################")
	print("t_0: {} Gyr".format(t0/1e9))
	print("e_0: {}".format(e0[ip_ind]))
	print("Omega_orb0: {} cyc/day".format(OmegaOrb0[ip_ind]))
	print("Omega_rot0: {} cyc/day".format(OmegaRot0[ip_ind]))

if __name__=="__main__":
	main()