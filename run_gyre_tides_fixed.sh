#!/bin/bash

# current parameters
m=$1
z=$2
job_n=${3-0} # take command line arg or else just label it as 0

# set up the directory for orbital evolution
# Check whether we have the base directory containing the photos
if [ ! -d "${cur_orbit_path}" ]; then
	# get MESA copy for filling in profiles
	cp -r ${cur_star_path} ${cur_orbit_path}
	# get glue scripts
	cp -r ${base_fidir}/orbev/* ${cur_orbit_path}
	# copy the initial orbital conditions and update-parameters
	cp ${base_fidir}/base_setup/gyre_base_setup/params* ${cur_orbit_path}
fi
cd ${cur_orbit_path}
rm ${cur_orbit_path}/photos/*
rm ${cur_orbit_path}/LOGS/history_full.data

# reformat the inlist to save every profile but no photos
sed -i "s/profile_interval=.*/profile_interval=1/g" inlist_MS # save every profile
sed -i "s/photo_interval=.*/photo_interval=100000/g" inlist_MS # photo interval is longer than max_model_number-> no photos
# save these profiles into the directory for this particular orbit
sed -i "s:log_directory=.*:log_directory='${orbev_fodir}':g" inlist_MS

# initial orbital configuration and simulation parameters
python params.py ${job_n}

# create a list of the original MESA photos so we can keep the photos directory clean
python create_photo_album.py ${cur_star_path} ${cur_orbit_path}

# if we aren't storing every profile, we need to generate the current chunk
# copy the correct photo or return an exit 1
# and edit the inlist_MS with the new max_model_num
python prepare_mesa_segment.py ${cur_star_path} ${cur_orbit_path} 0 ${job_n}
mesa_prep_exit_status=$?
if [ "${mesa_prep_exit_status}" -eq 0 ]; then
	echo "Restarting from current photo"
	# photo found successfully and inlist_MS modified
	# run the MESA model from current photo
	./re photo_cur
elif [ "${mesa_prep_exit_status}" -eq 1 ]; then
	echo "Restarting from ZAMS"
	# no preceding photo is available, so we just have to start from ZAMS .mod file
	./rn
elif [ "${mesa_prep_exit_status}" -eq 2 ]; then
	# no need to prepare any photos at all. We already have the profiles we need.
	echo "Profiles already exist. Continuing."
else
	echo "Not sure how you got here. Check your MESA segment preparation."
fi

# clean up the LOGS directory 
# don't confuse this for deleting profile*.data.GYRE, which we very much want to keep
rm ${orbev_fodir}/profile*.data
# consolidate profiles into a single file
python consolidate_profiles.py ${orbev_fodir}
# get rid of profile*.data.GYRE to prevent redundancy
rm ${orbev_fodir}/profile*.data.GYRE

# clean up photos directory
python clean_photo_album.py ${cur_star_path} ${cur_orbit_path}

# calculate the moments of inertia for the MESA models in the current chunk
python calculate_Is.py ${orbev_fodir}

# initialize orbital history file with initial orbital parameters
python initialize_state.py ${cur_orbit_path} ${cur_star} ${job_n}

# take some finite number of steps
for i in {1..10000}
do
	echo $i

	# generate new stellar profiles if needed via MESA simulation
	# check whether we are in the bounds of the currently generated profiles
	if ! python time_in_bounds.py ${orbev_fodir}; then
		echo "Running new segment of the MESA model at step $i"

		# clean up LOGS and photo directories
		rm ${orbev_fodir}/profile*
		rm ${orbev_fodir}/history.data

		python prepare_mesa_segment.py ${cur_star_path} ${cur_orbit_path} $i ${job_n}
		mesa_prep_exit_status=$?
		if [ "${mesa_prep_exit_status}" -eq 0 ]; then
			# photo found successfully and inlist_MS modified
			# run the MESA model from current photo
			./re photo_cur
		elif [ "${mesa_prep_exit_status}" -eq 1 ]; then
			# no preceding photo is available, so we just have to start from ZAMS .mod file
			./rn
		elif [ "${mesa_prep_exit_status}" -eq 2 ]; then
			# no need to prepare any photos at all. We already have the profiles we need.
			echo "Profiles already exist. Continuing."
		else
			echo "Not sure how you got here. Check your MESA segment preparation."
		fi

		# clean up the LOGS directory 
		# don't confuse this for deleting profile*.data.GYRE, which we very much want to keep
		rm ${orbev_fodir}/profile*.data
		# consolidate profiles into a single file
		python consolidate_profiles.py ${orbev_fodir}
		# get rid of profile*.data.GYRE to prevent redundancy
		rm ${orbev_fodir}/profile*.data.GYRE

		# clean up photos directory
		python clean_photo_album.py ${cur_star_path} ${cur_orbit_path}

		# calculate the moments of inertia for the MESA models in the current chunk
		python calculate_Is.py ${orbev_fodir}
	else
		echo "Continuing with step $i"
	fi

	# clean up the work directory
	if ((i>1)); then
		# get rid of the stellar profile before we get a fresh one
		# but only if we're making a fresh one!
		rm profile_cur.data.GYRE
	fi

	# take a timestep
	source ./GYRE_single_step.sh
done

# clean up output directory
rm ${orbev_fodir}/profile*
rm ${orbev_fodir}/history.data
rm ${orbev_fodir}/stellar_MOIs.txt

# get the history of orbital parameters
cp ${cur_orbit_path}/orbital_history.data ${orbev_fodir}/

# return to home dir
cd ${base_fidir}
