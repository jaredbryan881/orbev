import sys
import os

import numpy as np
from model_io import load_profile

# constants
G = 6.67430e-11  # [m3 kg-1 s-2]
Rsun = 695700000 # [m]
Msun = 1.989e30  # [kg]
au = 1.496e+11   # [m]

def main():
	base_finame=sys.argv[1]

	if os.path.exists("{}/stellar_MOIs.txt".format(base_finame)):
		print("Stellar MOI file already exists")
		return

	n_files=sum([len([file for file in files if ('.data.GYRE' in file)]) for root, dirs, files in sorted(os.walk(base_finame))])
	pinds = np.arange(1,1+n_files)
	Is=np.zeros(len(pinds))
	for (i,pind) in enumerate(pinds):
		profile, header = load_profile("{}/profile{}.data.GYRE".format(base_finame, pind))
		Is[i] = MOI(profile["M"], profile["r"]) # Msun*Rsun^2

	np.savetxt("{}/stellar_MOIs.txt".format(base_finame), Is)

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
