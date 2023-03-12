#!/bin/bash

# Base paths that won't change across stars/orbits

# set up MESA directories
export MESA_DIR=/home/jared/MIT/astero/mesa
export MESASDK_ROOT=/home/jared/MIT/astero/mesasdk
source $MESASDK_ROOT/bin/mesasdk_init.sh

# set up GYRE directories (with GYRE-tides enabled)
export GYRE_DIR=/home/jared/MIT/astero/gyre_v2/gyre

# path to directory containing orbev collection
base_fidir="/home/jared/MIT/astero/gyre_HATP2/orbev"

# path to directory we'll run the integration from
base_work_dir="/home/jared/MIT/astero/gyre_HATP2/orbev/work"

# path to the output directory
base_fodir="/home/jared/MIT/astero/gyre_HATP2/orbev/output"