import os
import numpy as np
from scipy.interpolate import interp1d

import ot

def select_profiles(cur_time, all_times, allowable_profiles, n_profiles):
	"""Get the indices of the N time samples in all_times closest to cur_time

	Arguments
	---------
	:param cur_time: float
		Current star age [yr]
	:param all_times: np.array
		Array of star ages from the MESA model [yr]
	:param allowable_profiles: list
		List of indices identifying the subset of all profiles which may be selected from
	:param n_profiles: int
		Number of nearest profiles to return

	Returns
	-------
	:return selected_pnums: list
		List of indices of nearest profiles in time to cur_time
	"""
	# get indices of each profile
	pnums=np.array([pnum for pnum in range(1,len(all_times)+1)])
	# mask out disallowed profiles
	allowable_mask=np.isin(pnums, allowable_profiles)

	# get closest n profiles in time
	dt_prof=np.abs(all_times - cur_time)
	# order indices by proximity to cur_time
	pnums_sort_inds=np.argsort(dt_prof)
	# filter sorted indices to just the allowable ones
	pnums_sort_inds=pnums_sort_inds[allowable_mask[pnums_sort_inds]]
	# take first n_profiles from the sorted indices
	selected_pnums=pnums_sort_inds[:n_profiles]
	# finally, just put those selected_pnums in order
	selected_pnums=np.sort(selected_pnums)
	
	return selected_pnums

def get_interpolation_axis(r1, r2):
	"""Create a new radial axis to allow interpolation between r1 and r2.
	Arguments
	---------
	:param r1: np.array
		Radial axis of the first stellar model
	:param r2: np.array
		Radial axis of the second stellar model
	Returns
	-------
	:return r1: np.array
		Radial axis of the first stellar model with new sampling
	:return r2: np.array
		Radial axis of the second stellar model with new sampling
	"""
	N1 = len(r1)
	N2 = len(r2)
	
	# Check whether resampling is needed
	if N1==N2:
		return r1, r2

	# Define uniform distribution over points in radial axis
	a=np.ones((N1,))/N1
	b=np.ones((N2,))/N2

	# TODO: combine these two cases by using dummy variables. Procedure is the same.
	
	# We want to resample both radii to the higher sampling case. So we check both cases.
	if N1>N2:
		# Calculate 1D transport map
		p = ot.emd_1d(r1, r2, a, b)
		# To prevent displacement of the stellar center and surface from their original locations
		# we round to two decimal places. This makes the map more 1:1 at the endpoints.
		p = np.around(N1*p, decimals=2)/N1

		# Create new radial axis for r2, where the position of the points in radius is given by
		# the average location of points transported from r1. I.e. barycentric projection
		r_interp=np.zeros(N1)
		for i in range(N1):
			nonzero = np.where(p[i,:]!=0)[0]
			for loc in nonzero:
				r_interp[i]+=N1*p[i,loc]*r2[loc]

		return r1, r_interp

	elif N2>N1:
		# Calculate 1D transport map
		p = ot.emd_1d(r2, r1, b, a)
		# To prevent displacement of the stellar center and surface from their original locations
		# we round to two decimal places. This makes the map more 1:1 at the endpoints.
		p=np.around(N2*p, decimals=2)/N2

		# Create new radial axis for r1, where the position of the points in radius is given by
		# the average location of points transported from r2. I.e. barycentric projection
		r_interp=np.zeros(N2)
		for i in range(N2):
			nonzero = np.where(p[i,:]!=0)[0]
			for loc in nonzero:
				r_interp[i]+=N2*p[i,loc]*r1[loc]

		return r_interp, r2

def lin_interp_2d(d1,d2,pct):
	"""Linearly interpolate some % between d1 and d2
	Arguments
	---------
	:param d1: np.array
		First data vector
	:param d2: np.array
		Second data vector
	:param pct: float
		Percentage of the way between d1 and d2
	Returns
	-------
	:return d_mid: np.array
		Data vector some % between d1 and d2
	"""
	return (1-pct)*d1 + (pct)*d2

def interpolate_single_quantity(r1, d1, r2, d2, pct, key):
	"""Interpolate between two curves with the same number of points.
	Arguments
	---------
	:param r: np.array
		Radius
	:param d1: np.array
		First data array
	:param d2: np.array
		Second data array
	:param pct: float
		A value between 0 and 1 representing the distance between d1 and d2
	Returns
	-------
	:return d_mid: np.array
		Data array pct between d1 and d2
	"""
	f=interp1d(d1["r"], d1[key])
	d1_interp=f(r1)
	f=interp1d(d2["r"], d2[key])
	d2_interp=f(r2)

	# simple linear interpolation
	d_mid=lin_interp_2d(d1_interp, d2_interp, pct)

	return d_mid