#!/bin/bash

# GYRE directory (with GYRE-tides enabled)
export GYRE_DIR=/home/jared/MIT/astero/gyre_v2/gyre

# path to directory containing orbev collection
base_fidir="/home/jared/MIT/astero/gyre_HATP2/orbev"

# path to directory we'll run code from
base_work_dir="/home/jared/MIT/astero/gyre_HATP2/orbev/work"

# path to place we'll save files to
base_fodir="/home/jared/MIT/astero/gyre_HATP2/orbev/output"

# current parameters
m=$1
z=$2
job_n=${3-0} # take command line arg or else just label it as 0

# create a place to work
# Path to MESA model for the current mass/metallicity
cur_star="M${m}_Z${z}"
cur_star_path="${base_work_dir}/${cur_star}"
# Path to MESA model for the current orbit, a copy of the one at cur_star_path
cur_orbit="orb${job_n}"
cur_orbit_path="${cur_star_path}_${cur_orbit}"
# Path to the MESA profiles and GYRE-tides output directory
fodir="${base_fodir}/${cur_star}"
orbev_fodir="${base_fodir}/${cur_star}_${cur_orbit}"
mkdir ${orbev_fodir}

# set up the directory for orbital evolution
# Check whether we have the base directory containing the photos
if [ ! -d "${cur_orbit_path}" ]; then
	# get MESA copy for filling in profiles
	cp -r ${cur_star_path} ${cur_orbit_path}
	# get glue scripts
	cp -r ${base_fidir}/orbev/* ${cur_orbit_path}
	# copy the initial orbital conditions and update-parameters
	cp ${base_fidir}/base_setup/gyre_base_setup/params.py ${cur_orbit_path}
fi
cd ${cur_orbit_path}

# reformat the inlist to save every profile but no photos
sed -i "s/profile_interval=.*/profile_interval=1/g" inlist_MS # save every profile
sed -i "s/photo_interval=.*/photo_interval=100000/g" inlist_MS # photo interval is longer than max_model_number-> no photos
# save these profiles into the directory for this particular orbit
sed -i "s:log_directory=.*:log_directory='${orbev_fodir}':g" inlist_MS

# info
python params.py ${job_n}

# create a list of the original MESA photos so we can keep the photos directory clean
python create_photo_album.py ${cur_orbit_path}

# if we aren't storing every profile, we need to generate the current chunk
# copy the correct photo or return an exit 1
# and edit the inlist_MS with the new max_model_num
python prepare_mesa_segment.py ${cur_orbit_path} 0
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
# clean up photos directory
python clean_photo_album.py ${cur_orbit_path}

# calculate the moments of inertia for the MESA models in the current chunk
python calculate_Is.py ${orbev_fodir}

# take some finite number of steps
for i in {1..10000}
do
	echo $i

	# get a fresh gyre inlist
	cp ${base_fidir}/base_setup/gyre_base_setup/gyre_tides.in .
	tide_foname="${orbev_fodir}/tide_orbit.h5"
	sed -i "s:summary_file\s=\s.*:summary_file = '${tide_foname}':g" gyre_tides.in

	# update the orbital parameters in the gyre inlist
	python update_orbital_parameters.py ${cur_orbit_path} $i $cur_star ${job_n} ${orbev_fodir}

	# generate new stellar profiles if needed via MESA simulation
	# check whether we are in the bounds of the currently generated profiles
	if ! python time_in_bounds.py ${orbev_fodir}; then
		echo "Running new segment of the MESA model at step $i"

		# clean up LOGS and photo directories
		rm ${orbev_fodir}/profile*
		rm ${orbev_fodir}/history.data
		rm ${cur_orbit_path}/photos/photo_cur

		python prepare_mesa_segment.py ${cur_orbit_path} $i
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
		rm ${orbev_fodir}/profile*.data
		# clean up photos directory
		python clean_photo_album.py ${cur_orbit_path}

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

	# update the stellar model in the working directory
	python update_stellar_profile.py ${orbev_fodir}

	# Calculate the orbtial evolution rates
	$GYRE_DIR/bin/gyre_tides gyre_tides.in

	# calculate secular rates of change and concatenate into a response history
	python calculate_orbev.py ${orbev_fodir} $i

	# break if no response was created
	repackage_status=$?
	if [ "${repackage_status}" -eq 0 ]; then
		echo "tide_orbit.h5 was created successfully at step $i"
	else
		echo "tide_orbit.h5 not created, exiting."
		break
	fi

	# clean up the current directory
	# get rid of the gyre inlist and tides_output.h5
	rm gyre_tides.in
	rm ${orbev_fodir}/tide_orbit.h5
done

# clean up output directory
rm ${orbev_fodir}/profile*
rm ${orbev_fodir}/history.data

cp ${cur_orbit_path}/orbital_history.data ${orbev_fodir}/

# return to home dir
cd ${base_fidir}
