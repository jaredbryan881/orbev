#!/bin/bash

source ./initialize_paths.sh

# current parameters
Ms=(1.36136 1.3634 1.3668 1.3702 1.3736 1.394 1.428 1.36 1.36 1.36 1.36 1.36 1.36 1.36)
Zs=(0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02262 0.02264262 0.02267655 0.0227331 0.02278965 0.0228462 0.0231855 0.023751)

job_n=${1-0} # take command line arg or else just label it as 0

m=${Ms[$job_n]}
z=${Zs[$job_n]}

photos=(x243000 x244000 x243000 x226000 x235000 x236000 x239000 x229000 x242000 x245000 x245000 x229000 x232000 x239000)
photo_str=${photos[$job_n]}

#############################
# Run background MESA model #
#############################
# path to the current MESA mode
cur_star="M${m}_Z${z}"
cur_star_path="${base_work_dir}/${cur_star}"

source ./continue_mesa_backbone.sh ${photo_str}
