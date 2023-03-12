#!/bin/bash

# current parameters
m=1.40
z=0.02
job_n=${1-0} # take command line arg or else just label it as 0

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
cur_orbit="orb${job_n}"
cur_orbit_path="${cur_star_path}_${cur_orbit}"
# Path to the MESA profiles and GYRE-tides output directory
orbev_fodir="${base_fodir}/${cur_star}_${cur_orbit}"
mkdir ${orbev_fodir}

# run the orbital evolution
source ./run_gyre_tides.sh $m $z $job_n
