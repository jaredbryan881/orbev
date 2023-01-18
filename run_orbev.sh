#!/bin/bash

Ms=(1.35)
Zs=(0.028)

for m in "${Ms[@]}"
do
	for z in "${Zs[@]}"
	do
		# run the backbone MESA model
		source ./run_mesa_backbone.sh $m $z

		# run the orbital evolution
		source ./run_gyre_tides.sh $m $z
	done
done