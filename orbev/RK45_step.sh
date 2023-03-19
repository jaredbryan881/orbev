#!/bin/bash

# take one timestep using the Runge-Kutta-Fehlberg method, aka RK4(5)
# involves six function evaluations
for rk_ind in {0..5}
do
	echo "RK4(5): function evaluation number ${rk_ind}"
	# get a fresh gyre inlist
	cp ${base_fidir}/base_setup/gyre_base_setup/gyre_tides.in .
	# Path to output file of gyre_tides
	tide_foname="${orbev_fodir}/tide_orbit_${rk_ind}.h5"
	# update location for output file
	sed -i "s:summary_file\s=\s.*:summary_file = '${tide_foname}':g" gyre_tides.in
	# update the orbital parameters in the inlist
	python update_inlist.py $i ${rk_ind} ${job_n} $m

	# Take a timestep using the Runge-Kutta-Fehlberg method, aka RK4(5)
	# update the stellar model in the working directory
	python update_stellar_profile.py ${orbev_fodir} ${rk_ind}

	# Calculate the orbtial evolution rates
	$GYRE_DIR/bin/gyre_tides gyre_tides.in

	# calculate secular rates of change and concatenate into a response history
	if ((rk_ind==0)); then
		# save the first one for long term storage
		python calculate_orbev.py ${orbev_fodir} tidal_response_history.h5 $rk_ind
	fi
	python calculate_orbev.py ${orbev_fodir} tidal_response_RKF_step.h5 $rk_ind

	# dimensionalize the orbital evolution rates
	# update the RKF_buffer data file with intermediate function evaluations
	python update_orbital_parameters.py ${cur_orbit_path} ${cur_star} ${job_n} ${orbev_fodir} ${rk_ind}

	# clean up the current directory
	# get rid of the gyre inlist and tides_output.h5
	rm gyre_tides.in
done

# clean up temporary files used to make the timestep
rm ${orbev_fodir}/tide_orbit*.h5
rm ${orbev_fodir}/tidal_response_RKF_step.h5
rm RKF_buffer.data