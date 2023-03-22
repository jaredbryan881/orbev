import numpy as np
import mesa_reader as mr
import sys

def main():
	cur_path=sys.argv[1]
	sh_finame="{}/LOGS/history_full.data"
	sh=mr.MesaData(sh_finame)

	np.savetxt("star_ages.txt", sh.star_age)

if __name__=="__main__":
	main()