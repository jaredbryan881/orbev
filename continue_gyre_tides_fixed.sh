#!/bin/bash

# current parameters
job_n=${1-0} # take command line arg or else just label it as 0

# Ms=(1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.48 1.24 1.36 1.36)
# Zs=(0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02996 0.01708)

Ms=(1.36 1.36136 1.3634 1.3668 1.3702 1.3736 1.394 1.428 1.36 1.36 1.36 1.36 1.36 1.36 1.36)
Zs=(0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02264262 0.02267655 0.0227331 0.02278965 0.0228462 0.0231855 0.023751)

select_n=$job_n%15

m=${Ms[$select_n]}
z=${Zs[$select_n]}

cur_star="M${m}_Z${z}"
cur_star_path="${base_work_dir}/${cur_star}"

cur_orbit="orb${job_n}_fixed_orbit"
cur_orbit_path="${cur_star_path}_${cur_orbit}"
# Path to the MESA profiles and GYRE-tides output directory
orbev_fodir="${base_fodir}/${cur_star}_${cur_orbit}"

cd ${cur_orbit_path}

# initial orbital configuration and simulation parameters
python params.py ${job_n}

# clean up LOGS and photo directories
rm ${orbev_fodir}/profile*
rm ${orbev_fodir}/history.data

nphotos=$(ls ${cur_star_path}/photos/ | wc -l)
if [ $nphotos -eq 0 ]
then
	echo "Photos somehow got deleted..."
	echo "Copying photos so we can continue"
	cp ${base_archive_dir}/work/${cur_star}/photos/* ${cur_star_path}/photos/
	nphotos=$(ls ${cur_star_path}/photos/ | wc -l)
	echo "Now have $nphotos photos. Continuing."
else
	echo "All $nphotos photos are here. Continuing."
fi

python prepare_mesa_segment.py ${cur_star_path} ${cur_orbit_path} 1 ${job_n}
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

#clean up photos directory
#python clean_photo_album.py ${cur_star_path} ${cur_orbit_path}

# calculate the moments of inertia for the MESA models in the current chunk
python calculate_Is.py ${orbev_fodir}

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

		nphotos=$(ls ${cur_star_path}/photos/ | wc -l)
		if [ $nphotos -eq 0 ]
		then
			echo "Photos somehow got deleted..." 
			echo "Copying photos so	we can continue"
			cp ${base_archive_dir}/work/${cur_star}/photos/* ${cur_star_path}/photos/
			nphotos=$(ls ${cur_star_path}/photos/ | wc -l)
			echo "Now have $nphotos	photos.	Continuing."
		else
			echo "All $nphotos photos are here. Continuing."
		fi

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
		#python clean_photo_album.py ${cur_star_path} ${cur_orbit_path}

		# calculate the moments of inertia for the MESA models in the current chunk
		python calculate_Is.py ${orbev_fodir}
	else
		echo "Continuing with step $i"
	fi

	rm profile_cur.data.GYRE

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
