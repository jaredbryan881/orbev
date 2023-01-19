#!/bin/bash

conda activate astero

# path to directory containing orbev collection
base_work_dir="/home/jared/MIT/astero/gyre_HATP2/orbev/"
# GYRE directory (with GYRE-tides enabled)
export GYRE_DIR=/home/jared/MIT/astero/gyre_v2/gyre

m=$1
z=$2
job_n=${3-0} # take command line arg or else just label it as 0

# create a place to work
cur_dir="M${m}_Z${z}"
cd work/${cur_dir}_${job_n}

# copy the update scripts
cp ../../orbev/*.py .

# copy the initial orbital conditions and update-parameters
cp ../../base_setup/gyre_base_setup/params.py .
python params.py $job_n

# create a list of the original photos so we can keep the photos directory clean
python create_photo_album.py

# if we aren't storing every profile, we need to generate the current chunk
# copy the correct photo or return an exit 1
# and edit the inlist_MS with the new max_model_num
if ! python prepare_mesa_segment.py 0; then
	# no preceding photo is available, so we just have to start from ZAMS .mod file
	./rn
else
	# run the MESA model from current photo
	./re photo_cur
fi
# clean up the LOGS directory
rm LOGS/profile*.data
# clean up photos directory
python clean_photo_album.py

# calculate the moments of inertia for the MESA models in the current chunk
python calculate_Is.py ${base_work_dir}/work/${cur_dir}_${job_n}/LOGS

# take some finite number of steps
for i in {1..10000}
do
	echo $i

	# get a fresh gyre inlist
	cp ${base_work_dir}/base_setup/gyre_base_setup/gyre_tides.in .

	if ((i>1)); then
		# calculate secular rates of change
		python calculate_orbev.py $i
	fi

	# update the orbital parameters in the gyre inlist
	python update_orbital_parameters.py $i $cur_dir $job_n

	# generate new stellar profiles if needed via MESA simulation
	# check whether we are in the bounds of the currently generated profiles
	if ! python time_in_bounds.py; then
		echo "Running new segment of the MESA model at step $i"

		# clean up LOGS and photo directories
		rm LOGS/profile*
		rm LOGS/history.data

		# copy the correct photo or return an exit 1
		# and edit the inlist_MS with the new max_model_num
		if ! python prepare_mesa_segment.py $i; then
			# no preceding photo is available, so we just have to start from ZAMS .mod file
			./rn
		else
			# run the MESA model from current photo
			./re photo_cur
		fi
		# clean up the LOGS directory
		rm LOGS/profile*.data
		# clean up photos directory
		python clean_photo_album.py

		# calculate the moments of inertia for the MESA models in the current chunk
		python calculate_Is.py ${base_work_dir}/work/${cur_dir}_${job_n}/LOGS
	else
		echo "Continuing with step $i"
	fi

	# clean up the work directory
	if ((i>1)); then
		# get rid of the stellar profile before we get a fresh one
		# but only if we're making a fresh one!
		rm profile_cur.data.GYRE
	fi

	# update the stellar model in the working directory
	python update_stellar_profile.py $i

	# Calculate the orbtial evolution rates
	$GYRE_DIR/bin/gyre_tides gyre_tides.in

	# clean up the current directory
	# get rid of the gyre inlist
	rm gyre_tides.in
	# move output files to profile directory
	mkdir profile$i
	mv tide_orbit.h5 profile$i/   # move forced oscillation summary

	# break if the directory is empty
	if [ "$(ls -A profile$i)" ]; then
		echo "tide_orbit.h5 was stored successfully for profile$i"
	else
		echo "tide_orbit.h5 not created, exiting."
		rmdir profile$i
		break
	fi
done

cd ..