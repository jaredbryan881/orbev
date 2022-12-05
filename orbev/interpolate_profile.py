import os
import numpy as np
from scipy import interpolate

import ot

def select_two_profiles(cur_time, all_times):
	"""Get the indices of t1 and t2 in all_times such that t1 <= cur_time < t2

	Arguments
	---------
	:param cur_time: float
		Current star age [yr]
	:param all_times: np.array
		Array of star ages from the MESA model [yr]
	"""
	pnum2 = np.argmin(all_times<cur_time)
	pnum1 = pnum2-1

	if cur_time<all_times[pnum1]:
		cur_time=all_times[pnum1]

	return pnum1, pnum2


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
	f=interpolate.interp1d(d1["r"], d1[key])
	d1_interp=f(r1)
	f=interpolate.interp1d(d2["r"], d2[key])
	d2_interp=f(r2)

	# simple linear interpolation
	d_mid=lin_interp_2d(d1_interp, d2_interp, pct)

	return d_mid
