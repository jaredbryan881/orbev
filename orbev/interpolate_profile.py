import os
import numpy as np
from scipy.interpolate import interp1d, make_interp_spline

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

def resample_radial_axes(rs):
	"""Create a new radial axis to allow interpolation between r1 and r2.

	The plan is to choose the profile with the most samples
	then approximate every other profile at the higher sampling
	by cascading Monge mappings (via barycentric projection of Kantorovich plans)

	Arguments
	---------
	:param rs: list
		Radial axis of each stellar model used in the current interpolation

	Returns
	-------
	:return rs_interp: list
		Radial axis with new sampling of each stellar model used in the current interpolation
	"""
	# Check whether resampling is needed
	Ns = [len(r) for r in rs]
	if len(set(Ns))==1:
		print("No resampling needed.")
		return rs
	print("Resampling {} profiles to {} samples".format(len(Ns), np.max(Ns)))

	# choose the profile with the most samples
	max_ind=np.argmax(Ns)

	# initialize list with the radial axis with the most samples
	rs_interp=[rs[max_ind]]

	# first cascade Monge maps backward in time
	for i in range(max_ind-1,-1,-1):
		r_monge = apply_transport_map(rs_interp[0], rs[i])
		rs_interp.insert(0,r_monge)

	# next cascade Monge maps forward in time
	for i in range(max_ind+1,len(Ns)):
		r_monge = apply_transport_map(rs_interp[-1], rs[i])
		rs_interp.append(r_monge)

	return rs_interp

def apply_transport_map(r1,r2):
	N1=len(r1)
	N2=len(r2)

	if N1==N2:
		return r2
	else:
		# define uniform distributions over radial samples
		a=np.ones((N1,))/N1
		b=np.ones((N2,))/N2

		# Calculate 1D transport map
		p = ot.emd_1d(r1, r2, a, b)

		# To prevent displacement of the stellar center and surface from their original locations
		# we round to two decimal places. This makes the map more 1:1 at the endpoints.
		p=np.around(N1*p, decimals=2)/N1

		# Create new radial axis for r2, where the position of the points in radius is given by
		# the average location of points transported from r1 to r2.
		# i.e. perform barycentric projection of the Kantorovich plan to approximate the Monge map 
		r_monge=np.zeros(N1)
		for i in range(N1):
			# consider contributions of points which accept nonzero mass
			for loc in np.where(p[i,:]!=0)[0]: 
				r_monge[i]+=N1*p[i,loc]*r2[loc]

		return r_monge

def resample_profiles(rs, ps):
	"""Resample radial profiles

	Arguments
	---------
	:param rs: list
		Radial axes to which the profiles will be resampled
	:param ps: list
		List of dicts, each dict has radial profiles of several quantities to be resampled

	Returns
	-------
	:return ps_interp: list
		List of dicts containing profiles now resampled to rs
	"""
	ps_interp=[]
	for (i,p) in enumerate(ps):
		p_interp_cur={}
		# resample each quantity's radial profile to the current r
		for key in p.keys():
			if key=="ind":
				p_interp_cur["ind"]=np.arange(1,len(rs[i])+1)
			elif key=="r":
				p_interp_cur["r"]=rs[i]
			else:
				f=interp1d(ps[i]["r"], ps[i][key])
				p_interp_cur[key]=f(rs[i])
		ps_interp.append(p_interp_cur)
	return ps_interp

def interp_2d(ps, cur_time, all_times, spline_order):
	"""Perform cubic spline interpolation on a 

	Arguments
	---------
	:param ps: list
		Data vectors at all_times
	:param cur_time: np.array
		First data vector
	:param all_times: np.array
		Second data vector
	:param spline_order: int
		B-spline degree

	Returns
	-------
	:return p_interp: np.array
		Data vector at cur_time
	"""
	# turn ps into a 2d array
	ps=np.array(ps)

	# perform cubic spline interpolation
	f = make_interp_spline(all_times, ps, k=spline_order, axis=0)

	# interpolate profile at cur_time
	p_interp=f(cur_time) 

	return p_interp