#!/bin/bash

source ./initialize_paths.sh

# current parameters
#Ms=(1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.48 1.24 1.36 1.36)
#Zs=(0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02996 0.01708)

Ms=(1.36 1.36136 1.3634 1.3668 1.3702 1.3736 1.394 1.428 1.36 1.36 1.36 1.36 1.36 1.36 1.36)
Zs=(0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02264262 0.02267655 0.0227331 0.02278965 0.0228462 0.0231855 0.023751)

job_n=${1-0} # take command line arg or else just label it as 0

select_n=$job_n%15
m=${Ms[$select_n]}
z=${Zs[$select_n]}

#############################
# Run background MESA model #
#############################
# path to the current MESA mode
cur_star="M${m}_Z${z}"
cur_star_path="${base_work_dir}/${cur_star}"

if [ ! -d ${cur_star_path} ]; then
	echo "Running MESA simulation."
	source ./run_mesa_backbone.sh $m $z
else
	echo "MESA model already exists."
fi

###########################
# Run orbital integration #
###########################
# Path to MESA model for the current orbit, a copy of the one at cur_star_path
cur_orbit="orb${job_n}_fixed_orbit"
cur_orbit_path="${cur_star_path}_${cur_orbit}"
# Path to the MESA profiles and GYRE-tides output directory
orbev_fodir="${base_fodir}/${cur_star}_${cur_orbit}"
mkdir ${orbev_fodir}

# run the orbital evolution
source ./run_gyre_tides_fixed.sh $m $z $job_n
