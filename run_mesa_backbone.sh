#!/bin/bash

# home directory
base_dir="/home/jared/MIT/astero/gyre_HATP2/orbev"

# current parameters
m=$1
z=$2

# create a place to work by copying the default setup
cur_dir="M${m}_Z${z}"
cp -r ${base_dir}/base_setup/mesa_base_setup ${base_dir}/work/${cur_dir}
cd ${base_dir}/work/${cur_dir}

# --- BEGIN PRE-MAIN SEQUENCE
echo "Running pre-main sequence evolution for a star of M=${m}Msun and Z=$z"
# get the correct initial mass and metallicity
sed -i "s/initial_mass=.*/initial_mass=$m/g" inlist_PMS
sed -i "s/initial_z=.*/initial_z=$z/g" inlist_PMS
sed -i "s/Zbase=.*/Zbase=$z/g" inlist_PMS

# run the pre-main sequence evolution to get a starting model at this mass/metallicity
./mk
./rn

# clean up the logs/photos directory, leaving only the .mod file
rm ${base_dir}/work/${cur_dir}/LOGS/*
rm ${base_dir}/work/${cur_dir}/photos/*
# --- END PRE-MAIN SEQUENCE

# --- BEGIN MAIN SEQUENCE
echo "Running main sequence evolution for a star of M=${m}Msun and Z=$z"
# get the correct initial mass and metallicity
sed -i "s/initial_mass=.*/initial_mass=$m/g" inlist_MS
sed -i "s/initial_z=.*/initial_z=$z/g" inlist_MS
sed -i "s/Zbase=.*/Zbase=$z/g" inlist_MS

# modify the inlist to set up the main-sequence evolution
sed -i 's/inlist_PMS/inlist_MS/g' inlist

# run the main-sequence with the sparese photo/profile sampling
./clean
./mk
./rn

# just get rid of the profiles to start fresh
rm ${base_dir}/work/${cur_dir}/LOGS/profile*

# rename the history file to protect it from being overwritten
mv ${base_dir}/work/${cur_dir}/LOGS/history.data ${base_dir}/work/${cur_dir}/LOGS/history_full.data
# --- END MAIN SEQUENCE

# reformat the inlist for subsequent runs
sed -i "s/profile_interval=.*/profile_interval=1/g" inlist_MS # save every profile
sed -i "s/photo_interval=.*/photo_interval=100000/g" inlist_MS # photo interval is longer than max_model_number-> no photos

# return to base directory
cd ${base_dir}
