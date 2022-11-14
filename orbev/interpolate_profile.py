import os
import numpy as np
from scipy import interpolate

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

def get_interpolation_axis(d1, d2, nbin=100):
	"""Create a new radius axis to allow interpolation between d1 and d2.

	Arguments
	---------
	:param d1: dict
		Dictionary of np.array data vectors
	:param d2: dict
		Dictionary of np.array data vectors
	:param nbin: int
		Number of bins over which the radius will be uniformly resampled

	Returns
	-------
	:return r1: np.array
		Radius of the first profile
	:return r2: np.array
		Radius of the second profile
	"""
	N1=len(d1["r"])
	N2=len(d2["r"])

	if N1!=N2:
		# we want to resample the curves so they have the same number of points
		# but the curves are variably sampled, so we want to respect sampling density.

		# First, let's calculate the density of the two curves
		bins=np.linspace(0,1,nbin+1)
		den1,_=np.histogram(d1["r"]/d1["r"].max(), bins)
		den2,_=np.histogram(d2["r"]/d2["r"].max(), bins)

		# Next, the easiest thing to do is just to resample both curves depending on whether
		# it is denser or less dense in this bin. This means that the interpolated curve will
		# always have more points than either of the end points.
		N_bin_rs = np.max([den1, den2], axis=0)

		r=np.zeros(N_bin_rs.sum())
		ind=0
		for i in range(nbin):
			r[ind:ind+N_bin_rs[i]]=np.linspace(bins[i], bins[i+1], N_bin_rs[i])
			ind+=N_bin_rs[i]
		# stretch these back to their respective lengths
		r1=r*d1["r"].max()
		r2=r*d2["r"].max()

	else:
		# no interpolation necessary
		r1=d1["r"]
		r2=d2["r"]

	return r1, r2

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
	d_mid=np.zeros(len(r1))

	f=interpolate.interp1d(d1["r"], d1[key])
	d1_interp=f(r1)
	f=interpolate.interp1d(d2["r"], d2[key])
	d2_interp=f(r2)

	# simple linear interpolation
	d_mid=lin_interp_2d(d1_interp, d2_interp, pct)

	return d_mid
