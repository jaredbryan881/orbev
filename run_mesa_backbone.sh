#!/bin/bash

# home directory
orbev_dir="/home/jared/MIT/astero/gyre_HATP2/orbev"
# directory containing the MESA stellar profiles
base_mesa_work_dir="/home/jared/MIT/astero/mesa_HATP2/live_planet"
# GYRE directory (with GYRE-tides enabled)
export MESA_DIR=/home/jared/MIT/astero/mesa

m=$1
z=$2
job_n=${3-0}

# create a place to work by copying the default setup
cur_dir="M${m}_Z${z}"
cp -r ${base_mesa_work_dir}/base_setup output/${cur_dir}_${job_n}
cd output/${cur_dir}_${job_n}

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
rm LOGS/*
rm photos/*

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
rm LOGS/profile*

# rename the history file to protect it from being overwritten
mv LOGS/history.data LOGS/history_full.data

# --- END MAIN SEQUENCE

# return to base directory
cd ../..

# link photos to their indices in history.data
python link_photos_to_history.py