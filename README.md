This package provides tools for coupling the MESA stellar evolution code and the GYRE-tides stellar pulsation code.
The process of coupled stellar and orbital evolution is as follows:
1. Run the MESA stellar evolution code with a short timestep and sparse photo and profile saving interval.
2. Choose an initial orbital configuration and stellar age.
3. Restart MESA from a photo near the chosen stellar age and save profiles at short sampling intervals over a given time period.
4. Run GYRE-tides to calculate the orbital evolution rates, update the orbital elements, and update the stellar model.
5. Repeat (4) until the system age evolves past the range of currently stored profiles.
6. Run (3) to create a new chunk of saved profiles and repeat (5).

Integration of the orbital evolution rates is performed with the Runge-Kutta-Fehlberg method.
This is a fully explicit 4th order method with 5th order error estimation, which we use to adaptively set the timestep.
Evaluation of the orbital evolution rates at arbitrary times is supported, with optimal transport-based interpolation between stellar models enabling the MESA and GYRE-tides simulations to be decoupled.
You can also perform both forward and time-reversed simulations.

Plotting scripts display the simulation outputs, including snapshots of the evolving stellar structure, tidal migration trajectories, and time-mode decompositions of the photometric amplitudes and orbital evolution rates.
Abbreviated simulation outputs for reproducing the figures in Bryan et al. (2024) can be found on Zenodo (10.5281/zenodo.11509036).
