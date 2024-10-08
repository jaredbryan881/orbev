import glob
import sys
import os

def main():
	# get base location of MESA photos
	cur_star_path=sys.argv[1]
	cur_orbit_path=sys.argv[2]

	# read photo names
	photo_list = glob.glob("{}/photos/*".format(cur_star_path))

	# see whether photo album already exists
	if os.path.isfile("{}/photo_album.txt".format(cur_orbit_path)):
		return

	# write photo names to file for safe keeping
	with open("{}/photo_album.txt".format(cur_orbit_path), "w") as f:
		for photo_name in photo_list:
			f.write('%s\n' % photo_name)
		print("Photo album has been created")

if __name__=="__main__":
	main()
