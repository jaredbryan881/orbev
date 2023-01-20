#!/bin/bash

conda activate astero

# GYRE directory (with GYRE-tides enabled)
export GYRE_DIR=/home/jared/MIT/astero/gyre_v2/gyre

# path to directory containing orbev collection
base_dir="/home/jared/MIT/astero/gyre_HATP2/orbev"

# current parameters
m=$1
z=$2
job_n=${3-0} # take command line arg or else just label it as 0

# create a place to work
cur_dir="M${m}_Z${z}"
cur_path="${base_dir}/work/${cur_dir}"
cd work/${cur_dir}

# create a place to do GYRE-tides work
cp -r ${base_dir}/orbev ${cur_dir}_${job_n}
cd ${cur_dir}_${job_n}

# copy the initial orbital conditions and update-parameters
cp ${base_dir}/base_setup/gyre_base_setup/params.py .
python params.py ${job_n}

# create a list of the original MESA photos so we can keep the photos directory clean
python create_photo_album.py ${cur_path}

# if we aren't storing every profile, we need to generate the current chunk
# copy the correct photo or return an exit 1
# and edit the inlist_MS with the new max_model_num
python prepare_mesa_segment.py ${cur_path} 0
mesa_prep_exit_status=$?
if [ "${mesa_prep_exit_status}" -eq 0 ]; then
	echo "Restarting from current photo"
	# photo found successfully and inlist_MS modified
	# run the MESA model from current photo
	cd ${cur_path}
	./re photo_cur
	cd ${cur_path}/${cur_dir}_${job_n}
elif [ "${mesa_prep_exit_status}" -eq 1 ]; then
	echo "Restarting from ZAMS"
	# no preceding photo is available, so we just have to start from ZAMS .mod file
	cd ${cur_path}
	./rn
	cd ${cur_path}/${cur_dir}_${job_n}
elif [ "${mesa_prep_exit_status}" -eq 2 ]; then
	# no need to prepare any photos at all. We already have the profiles we need.
	echo "Profiles already exist. Continuing."
else
	echo "Not sure how you got here. Check your MESA segment preparation."
fi

# clean up the LOGS directory 
# don't confuse this for deleting profile*.data.GYRE, which we very much want to keep
rm ${cur_path}/LOGS/profile*.data
# clean up photos directory
python clean_photo_album.py ${cur_path}

# calculate the moments of inertia for the MESA models in the current chunk
python calculate_Is.py ${cur_path}/LOGS

# take some finite number of steps
for i in {1..10000}
do
	echo $i

	# get a fresh gyre inlist
	cp ${base_dir}/base_setup/gyre_base_setup/gyre_tides.in .

	if ((i>1)); then
		# calculate secular rates of change
		python calculate_orbev.py $i
	fi

	# update the orbital parameters in the gyre inlist
	python update_orbital_parameters.py ${cur_path} $i $cur_dir ${job_n}

	# generate new stellar profiles if needed via MESA simulation
	# check whether we are in the bounds of the currently generated profiles
	if ! python time_in_bounds.py ${cur_path}; then
		echo "Running new segment of the MESA model at step $i"

		# clean up LOGS and photo directories
		rm ${cur_path}/LOGS/profile*
		rm ${cur_path}/LOGS/history.data

		python prepare_mesa_segment.py ${cur_path} $i
		mesa_prep_exit_status=$?
		if [ "${mesa_prep_exit_status}" -eq 0 ]; then
			# photo found successfully and inlist_MS modified
			# run the MESA model from current photo
			cd ${cur_path}
			./re photo_cur
			cd ${cur_path}/${cur_dir}_${job_n}
		elif [ "${mesa_prep_exit_status}" -eq 1 ]; then
			# no preceding photo is available, so we just have to start from ZAMS .mod file
			cd ${cur_path}
			./rn
			cd ${cur_path}/${cur_dir}_${job_n}
		elif [ "${mesa_prep_exit_status}" -eq 2 ]; then
			# no need to prepare any photos at all. We already have the profiles we need.
			echo "Profiles already exist. Continuing."
		else
			echo "Not sure how you got here. Check your MESA segment preparation."
		fi

		# clean up the LOGS directory
		rm ${cur_path}/LOGS/profile*.data
		# clean up photos directory
		python clean_photo_album.py ${cur_path}

		# calculate the moments of inertia for the MESA models in the current chunk
		python calculate_Is.py ${cur_path}/LOGS
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
	python update_stellar_profile.py ${cur_path}

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

cd ${base_dir}