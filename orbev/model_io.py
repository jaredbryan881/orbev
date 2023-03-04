import numpy as np
import pandas as pd

import os
import re
import fileinput

def save_profile(profile, header):
	"""Save interpolated profile.

	Argumnets
	---------
	:param profile: dict
		Dictionary of radial profiles of parameters
	:param header: list
		Metadata for the profile
	"""
	ncol=len(profile.keys())
	npts=len(profile["ind"])
	profile_data=np.zeros((npts, ncol))
	for (i,key) in enumerate(profile.keys()):
		profile_data[:,i]=profile[key]

	with open("profile_cur.data.GYRE", "a+") as f:
		# formatting is to match profile.data.GYRE
		# https://gyre.readthedocs.io/en/stable/ref-guide/stellar-models/mesa-file-format-v1.01.html
		# header has two leading spaces, an integer number of points, 16-decimal exponential notation floats,
		# and another integer for the version number, each separated by 5 spaces
		# body has 5 leading spaces, an integer index, and 16-decimal place exponential floats,
		# each separated by 5 spaces
		np.savetxt(f, header, fmt='  %i     %1.16e     %1.16e     %1.16e     %i')
		np.savetxt(f, profile_data, fmt='     %i' + '     %1.16e'*(ncol-1))

def load_profile(finame):
	"""Load a profile.data.GYRE file

	Arguments
	---------
	:param finame: str
		Path to the profile.data.GYRE file

	Returns
	-------
	:return profile: dict
		Dictionary of np.array data vectors
	:return header:
	"""
	header=np.loadtxt(finame, max_rows=1) # Npts, Mass [g], Radius [cm], Luminosity [erg/s]
	data=np.loadtxt(finame, skiprows=1, max_rows=int(header[0]))
	profile={"ind"       : data[:,0],  # grid point index
			 "r"         : data[:,1],  # radial coordinate [cm]
			 "M"         : data[:,2],  # interior mass [g]
			 "L"         : data[:,3],  # interior luminosity [erg/s]
			 "P"         : data[:,4],  # total pressure [Ba]
			 "T"         : data[:,5],  # temperature [K]
			 "rho"       : data[:,6],  # density [g/cm^3]
			 "nab"       : data[:,7],  # temperature gradient
			 "N2"        : data[:,8],  # Brunt-Vaisala Frequency [1/s^2]
			 "Gam1"      : data[:,9],  # adiabatic exponent
			 "nab_ad"    : data[:,10], # adiabatic temmperature gradient
			 "nu_T"      : data[:,11], # thermodynamic coefficient
			 "kap"       : data[:,12], # opacity
			 "kap_T"     : data[:,13], # opacity partial
			 "kap_rho"   : data[:,14], # opacity partial
			 "eps_n"     : data[:,15], # nuclear energy generation rate
			 "eps_n_T"   : data[:,16], # nuclear energy generation rate partial
			 "eps_n_rho" : data[:,17], # nuclear energy generation rate partial
			 "Omega_rot" : data[:,18]  # Rotational velocity
			}

	return profile, header

def load_stellar_state(finame):
	"""Read the current mass and radius of the star.

	Arguments
	---------
	:param finame: str
		Input filename for current profile

	Returns
	-------
	:return cur_R: float
		current radius
	:return cur_M: float
		current mass
	"""
	data=np.loadtxt(finame, max_rows=1)
	cur_M=data[1]
	cur_R=data[2]

	return cur_R, cur_M

def load_orbital_state(oh_finame):
	"""Read the current orbital configuration

	Arguments
	---------
	:param oh_finame: str
		Orbital history file name

	Returns
	-------
	:return time: float
		Current stellar age
	:return a: float
		Current orbital semi-major axis [au]
	:return e: float
		Current orbital eccentricity []
	:return OmegaRot: float
		Current spin frequency [cyc/day]
	"""
	data=pd.read_csv(oh_finame)

	return data["time"].values[-1], data["a"].values[-1], data["e"].values[-1], data["OmegaRot"].values[-1]

def update_history(time, a, e, OmegaRot):
	"""Update orbital history file with current orbital configuration

	Arguments
	---------
	:param time: float
		Current stellar age
	:param a: float
		Current orbital semi-major axis [au]
	:param e: float
		Current orbital eccentricity []
	:param OmegaRot: float
		Current spin frequency [cyc/day]
	"""
	# check if it exists yet
	if not os.path.exists("orbital_history.data"):
		# initialize file with headers
		header="time,a,e,OmegaRot"
	else:
		header=''

	# append configuration
	with open("orbital_history.data", "a+") as f:
		np.savetxt(f, np.array([[time,a,e,OmegaRot]]), delimiter=',', header=header, comments='')

def update_orbital_parameters(OmegaOrb, OmegaRot, e, finame="gyre_tides.in"):
	"""Update orbital parameters in the GYRE-tides input file.

	Arguments
	---------
	:param OmegaOrb: float
		New orbital frequency [cyc/day]
	:param OmegaRot: float
		New rotational frequency [cyc/day]
	:param e: float
		New orbital eccentricity
	"""
	with open(finame, "r+") as file:
		text = file.read()
		text = re.sub(r"Omega_orb=.*", "Omega_orb={}".format(OmegaOrb), text)
		text = re.sub(r"Omega_rot=.*", "Omega_rot={}".format(OmegaRot), text)
		text = re.sub(r"e=.*", "e={}".format(e), text)

		file.seek(0)
		file.write(text)
		file.truncate()
