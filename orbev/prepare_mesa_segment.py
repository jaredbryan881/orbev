import numpy as np
import mesa_reader as mr

import re
import sys
import glob
import shutil

import params

from model_io import load_orbital_state

def main():
	pind=int(sys.argv[1])
	# get the current time
	if pind==0:
		# load initial time
		cur_time = params.t0
	else:
		# load current time
		oh_finame="orbital_history.data"
		cur_time,_,_,_=load_orbital_state(oh_finame)

	# locate the nearest photo
	photo_string, photo_num = get_nearest_photo(cur_time)

	# TODO: deal with earliest segment (what do I do when there is no photo?)
	# copy the photo
	shutil.copy("./photos/{}".format(photo_string), "./photos/photo_cur")

	# update the MESA inlist's max_model_num
	with open("inlist_MS", "r+") as file:
		text = file.read()
		text = re.sub(r"max_model_number=.*", "max_model_number={}".format(photo_num+105), text)

		file.seek(0)
		file.write(text)
		file.truncate()

def get_nearest_photo(cur_time):
	"""Find the nearest preceding photo to cur_time

	Arguments
	---------
	:param cur_time: float
		Current stellar age [yr]

	Returns
	-------
	:return photo_string: str
		Name of the nearest preceding photo
	:return photo_num: int
		Model number of the nearest preceding photo
	"""
	# Read stellar history file
	sh_finame="./LOGS/history_full.data"
	sh=mr.MesaData(sh_finame)

	# get a list of profiles
	photo_list = glob.glob("./photos/*")
	# keep just the photo names without the full path
	photo_list = [photo.split('/')[-1] for photo in photo_list]
	# interpret the photo name as an integer for association with model_number
	photo_model_num = [int(photo.split('x')[-1]) for photo in photo_list]
	photo_model_ind = [pmn-sh.model_number[0] for pmn in photo_model_num]

	# we want to get the closest photo to cur_time, but also earlier than cur_time
	# first let's get the stellar age at each photo
	photo_ages = np.array([sh.star_age[photo_num-1] for photo_num in photo_model_ind])
	# then let's get the time difference between each photo and cur_time
	photo_dt = cur_time - photo_ages
	# check if we are in the first segment, which has no preceding photo
	if all(photo_dt<0):
		return None
	# otherwise select the photo index and corresponding photo string
	photo_ind = min((p,i) for i,p in enumerate(photo_dt) if p>0)[1]
	photo_string = photo_list[photo_ind]

	return photo_string, photo_model_num[photo_ind]

if __name__=="__main__":
	main()