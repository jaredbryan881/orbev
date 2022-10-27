#!/bin/bash

# directory containing template gyre inlists
base_dir="/home/jared/MIT/astero/gyre_HATP2/parameter_search/base_setup"
# directory containing the MESA stellar profiles
mesa_dir="/home/jared/MIT/astero/mesa_HATP2/live_planet"
# GYRE directory (with GYRE-tides enabled)
export GYRE_DIR=/home/jared/MIT/astero/gyre_v2/gyre

Zs=(0.02)
Ms=(1.40)

for m in "${Ms[@]}"
do
	for z in "${Zs[@]}"
	do
		# create a place to work
		cur_dir="M${m}_Z${z}"
		mkdir $cur_dir
		cd $cur_dir

		# copy the update scripts
		cp ../update_orbital_parameters.py .
		cp ../update_stellar_profile.py .
		cp ../interpolate_profile.py .
		cp ../calculate_orbev.py .
		cp ../unit_conversion.py .
		cp ../calculate_Is.py .
		cp ../model_io.py .

		# copy the parameters
		cp ../params.py .

		# copy the first profile into working directory as our starting point
		cp $mesa_dir/$cur_dir/LOGS/profile1.data.GYRE profile_cur.data.GYRE

		# calculate the moments of inertia for the (noninterpolated MESA models)
		python calculate_Is.py $cur_dir

		# take some finite number of steps
		for i in {1..200}
		do
			echo $i

			# get a fresh gyre inlist
			cp $base_dir/gyre_tides.in .

			if ((i>1)); then
				# calculate secular rates of change
				python calculate_orbev.py $i $cur_dir
			fi

			# update the orbital parameters in the gyre inlist
			python update_orbital_parameters.py $i $cur_dir

			if ((i>1)); then
				# get rid of the stellar profile before we get a fresh one
				# but only if we're making a fresh one!
				rm profile_cur.data.GYRE
			fi
			# update the stellar model in the working directory
			python update_stellar_profile.py $i $cur_dir

			# Calculate the orbtial evolution rates
			$GYRE_DIR/bin/gyre_tides gyre_tides.in

			# clean up the current directory
			# get rid of the gyre inlist
			rm gyre_tides.in

			mkdir profile$i
			mv *.h5 profile$i/
		done

		cd ..
	done
done
