import numpy as np
import glob

def main():
	# read photo names
	photo_list = glob.glob("./photos/*")
	# write their names to file for safe keeping
	with open("photo_album.txt", "w") as f:
		for photo_name in photo_list:
			f.write('%s\n' % photo_name)
		print("Photo album has been created")

if __name__=="__main__":
	main()