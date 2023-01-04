This package provides tools for coupling the MESA stellar evolution code and the GYRE-tides stellar pulsation code.
First, you run the MESA stellar evolution code to obtain a list of stellar profiles.
You then choose an initial orbital configuration and initial stellar model.
You then run GYRE-tides to calculate the orbital evolution rates.
Arbitrary timesteps are supported, with optimal transport-based interpolation between stellar models enabling the MESA and GYRE-tides simulations to be decoupled.
You can also perform both forward and time-reversed simulations.