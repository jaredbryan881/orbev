#!/bin/bash

# current parameters
m=$1
z=$2

# create a place to work by copying the default setup
cp -r ${base_fidir}/base_setup/mesa_base_setup ${cur_star_path}
cd ${cur_star_path}

# --- BEGIN PRE-MAIN SEQUENCE
echo "Running pre-main sequence evolution for a star of M=${m}Msun and Z=$z"
# set the initial mass and metallicity
sed -i "s/initial_mass=.*/initial_mass=$m/g" inlist_PMS
sed -i "s/initial_z=.*/initial_z=$z/g" inlist_PMS
sed -i "s/Zbase=.*/Zbase=$z/g" inlist_PMS

# run the pre-main sequence evolution to get a starting model at this mass/metallicity
./mk
./rn

# clean up the logs/photos directory, leaving only the .mod file
rm ${cur_star_path}/LOGS/*
rm ${cur_star_path}/photos/*
# --- END PRE-MAIN SEQUENCE

# --- BEGIN MAIN SEQUENCE
echo "Running main sequence evolution for a star of M=${m}Msun and Z=$z"
# set the initial mass and metallicity
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
rm ${cur_star_path}/LOGS/profile*

# rename the history file to protect it from being overwritten
mv ${cur_star_path}/LOGS/history.data ${cur_star_path}/LOGS/history_full.data
# extract the times so we don't have to store the full history file for each orbital simulation
python extract_times.py ${cur_star_path}
# --- END MAIN SEQUENCE

# return to base directory
cd ${base_fidir}
