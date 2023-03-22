import numpy as np
import mesa_reader as mr
import sys

def main():
	cur_path=sys.argv[1]
	sh_finame="{}/LOGS/history_full.data".format(cur_path)
	sh=mr.MesaData(sh_finame)

	# save model numbers and associated stellar ages
	np.savetxt("{}/LOGS/star_ages.txt".format(cur_path), [sh.model_number, sh.star_age])

if __name__=="__main__":
	main()