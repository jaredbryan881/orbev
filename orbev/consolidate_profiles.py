import numpy as np
import sys

from model_io import load_profile
import pickle as pkl

import params

def main():
	base_profile_dir=sys.argv[1]

	headers  = []
	profiles = []
	for pnum in params.allowable_profiles:
		try:
			profile,header=load_profile("{}/profile{}.data.GYRE".format(base_profile_dir, pnum))
		except:
			continue
		headers.append(header)
		profiles.append(profile)

	with open("{}/profiles.pkl".format(base_profile_dir), "wb") as f:
		pkl.dump([headers, profiles], f)

if __name__=="__main__":
	main()