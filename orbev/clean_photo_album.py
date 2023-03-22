import os
import glob
import numpy as np

import sys

def main():
	cur_star_path=sys.argv[1]
	cur_orbit_path=sys.argv[2]
	# read current photo names
	photo_list_all = glob.glob("{}/photos/*".format(cur_star_path))
	# read photo names we'd like to keep
	with open("{}/photo_album.txt".format(cur_orbit_path), "r") as f:
		photo_list_keep=f.read().split('\n')

	# remove the photos not in photo_list_keep
	for photo in photo_list_all:
		if photo not in photo_list_keep:
			os.remove(photo)

if __name__=="__main__":
	main()
