#!/bin/bash

export MESA_DIR=/home/jared/MIT/astero/mesa
export MESASDK_ROOT=/home/jared/MIT/astero/mesasdk
source $MESASDK_ROOT/bin/mesasdk_init.sh
export GYRE_DIR=/home/jared/MIT/astero/gyre_v2/gyre


Ms=(1.35 1.35)
Zs=(0.028 0.028)

job_n=${1-0}

base_work_dir="."

for m in "${Ms[@]}"
do
	for z in "${Zs[@]}"
	do
		cur_star=M${m}_Z${z}

		if [ ! -d "${base_work_dir}/work/${cur_star}" ]; then
			echo "Running MESA simulation."
			# run the backbone MESA model
			source ./run_mesa_backbone.sh $m $z
		else
			echo "MESA model already exists."
		fi

		# run the orbital evolution
		source ./run_gyre_tides.sh $m $z $job_n
	done
done
