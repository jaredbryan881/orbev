import sys
import os

import pickle as pkl

import numpy as np

# constants
G = 6.67430e-11  # [m3 kg-1 s-2]
Rsun = 695700000 # [m]
Msun = 1.989e30  # [kg]
au = 1.496e+11   # [m]

def main():
	base_profile_dir=sys.argv[1]

	print("Calculating Stellar MOIs")
	with open("{}/profiles.pkl".format(base_profile_dir), "rb") as f:
		headers,profiles=pkl.load(f)

	Is=np.zeros(len(profiles))
	for (i,profile) in enumerate(profiles):
		Is[i] = MOI(profile["M"], profile["r"]) # Msun*Rsun^2

	np.savetxt("{}/stellar_MOIs.txt".format(base_profile_dir), Is)

def MOI(M, R):
	"""Calculate the moment of inertia of a 1D stellar model as the sum of
	moments of inertia in a series of thick spherical shells.

	Arguments
	---------
	:param M: np.array
		Mass a function of radius [g]
	:param R: np.array
		Radius [cm]

	Returns
	-------
	:return I: float
		Moment of inertia [Msun * Rsun**2]
	"""
	M/=1000 # [kg]
	M/=Msun # [Msun]

	R/=100  # [m]
	R/=Rsun # [Rsun]

	# mass in each layer
	dM = M[1:]-M[:-1]

	# moment of inertia in each layer
	I = (2/5)*dM*((R[1:]**5 - R[:-1]**5)/(R[1:]**3 - R[:-1]**3))

	# total moment of inertia
	I = np.sum(I)

	return I

if __name__=="__main__":
	main()
