#!/bin/bash

# current parameters
photo_str=$1

# create a place to work by copying the default setup
cd ${cur_star_path}

# --- CONTINUE MAIN SEQUENCE
echo "Running main sequence evolution for a star of M=${m}Msun and Z=$z"
./re ${photo_str}

# just get rid of the profiles to start fresh
rm ${cur_star_path}/LOGS/profile*

# rename the history file to protect it from being overwritten
mv ${cur_star_path}/LOGS/history.data ${cur_star_path}/LOGS/history_full.data
# extract the times so we don't have to store the full history file for each orbital simulation
python extract_times.py ${cur_star_path}
# --- END MAIN SEQUENCE

# return to base directory
cd ${base_fidir}
