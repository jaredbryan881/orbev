#!/bin/bash

Ms=(1.35 1.35)
Zs=(0.028 0.028)

job_n=${1-0}

for m in "${Ms[@]}"
do
	for z in "${Zs[@]}"
	do
		cur_star=M${m}_Z${z}

		if [ ! -d "./work/${cur_star}" ]; then
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