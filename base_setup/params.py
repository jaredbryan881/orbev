import numpy as np

###########
# Toggles #
###########
live_orbit=True

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
mesa_diname="/home/jared/MIT/astero/gyre_HATP2/benchmarking/profiles"
