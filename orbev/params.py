#################################
# Initial Orbital Configuration #
#################################
e0 = 0.5        # []
OmegaOrb0 = 0.177 # [cyc/day]
OmegaRot0 = OmegaOrb0*4 # [cyc/day]

#################
# Update limits #
#################
max_dt=1e5         # [yr]
max_de=1e-5        # []
max_da=1e-5        # [au]
max_dOmegaRot=1e-5 # [cyc/day]

###################
# GYRE Param File #
###################
#finame="./gyre_orbit.in"
finame="./gyre_tides.in"