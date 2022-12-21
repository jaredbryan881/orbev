import numpy as np

###########
# Toggles #
###########
live_orbit=False

#################################
# Initial Orbital Configuration #
#################################
e0 = np.array([0.5])
OmegaOrb0 = np.array([0.1774])
OmegaRot0 = np.array([0.1774])*4.8

#################
# Update limits #
#################
max_dt=1e5         # [yr]
max_de=1e-5        # []
max_da=1e-5        # [au]
max_dOmegaRot=1e-5 # [cyc/day]

#######################
# Paths and Filenames #
#######################
gyre_inlist="./gyre_tides.in"
#gyre_inlist="./gyre_orbit.in"
mesa_diname="/home/jared/MIT/astero/gyre_HATP2/benchmarking/profiles"

# subset of profiles which can be selected from
allowable_profiles=np.arange(11700,11800,2)

def main():
	import sys
	ip_ind=int(sys.argv[1])
	print("####################################")
	print("#### Initial Orbital Parameters ####")
	print("####################################")
	print("e_0: {}".format(e0[ip_ind]))
	print("Omega_orb0: {} cyc/day".format(OmegaOrb0[ip_ind]))
	print("Omega_rot0: {} cyc/day".format(OmegaRot0[ip_ind]))

if __name__=="__main__":
	main()