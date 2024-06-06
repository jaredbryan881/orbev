import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import mesa_reader as mr
import os

import warnings
warnings.filterwarnings("ignore")

Rsun = 695700000 # [m]
au = 1.496e+11   # [m]

def main():
	# plot subsynchronous data
	suffix="loworot"
	plot_orbev_equilibrium(suffix)

	# plot supersynchronous data
	suffix="highorot"
	plot_orbev_equilibrium(suffix)

def plot_orbev_equilibrium(suffix):
	base_diname="/home/jaredbryan/MIT/astero/gyre_HATP2/orbev/data/long_time/"
	base_star=["M1.36_Z0.02262", "M1.48_Z0.02262", "M1.24_Z0.02262", "M1.36_Z0.02996", "M1.36_Z0.01708"]

	colors=['k'] + 3*['green'] + 3*['steelblue'] + 3*['crimson'] + 2*['purple'] + 2*['chocolate'] + 2*["slategray"]
	linestyles=['solid'] + 3*['solid', 'dashed', (0, (1,0.5))] + 3*['solid', "dashed"]
	linewidths=16*[2]

	orbev_finame="tidal_response_history.h5"
	history_finame="orbital_history.data"

	# Initialize plot
	fig,axs=plt.subplots(6,2,sharex=True,sharey='col',figsize=(7,10))
	axs[-1,0].set_xlabel("Time [Myr]", fontsize=14)
	axs[-1,1].set_xlabel("Time [Myr]", fontsize=14)

	axs[0,0].set_xlim(-1.5, 1.7)

	#axs[0,0].set_ylim(0.3, 0.62)
	#axs[0,1].set_ylim(0.055, 0.09)

	# base trajectory
	time_base, e_base, a_base, e_dot_freq_base, a_dot_freq_base = process_tracks(base_diname, base_star[0], 0, suffix, history_finame, orbev_finame)
	t0=2.6e9
	t0s=[2.6e9]*10 + [3.1e9, 2.1e9] + [2.6e9]*4
	t_ind=np.argmin(np.abs(time_base-t0))
	e0=e_base[t_ind]
	a0=a_base[t_ind]
	t0=time_base[t_ind]

	import copy
	a_set=copy.deepcopy(a0)
	e_eq_base, a_eq_base = equilibrium_evolution(time_base, e_dot_freq_base, a_dot_freq_base, e0=e0, a0=a0, t0=t0)
	add_trajectory(axs[0], time_base-t0, e_base, a_base, a_set=a_set)
	add_trajectory(axs[1], time_base-t0, e_base, a_base, a_set=a_set)
	add_trajectory(axs[2], time_base-t0, e_base, a_base, a_set=a_set)
	add_trajectory(axs[3], time_base-t0, e_base, a_base, a_set=a_set)
	add_trajectory(axs[4], time_base-t0, e_base, a_base, a_set=a_set)
	add_trajectory(axs[5], time_base-t0, e_base, a_base, a_set=a_set)

	# variation in eccentricity
	for i in range(1,4):
		time, e, a, e_dot_freq, a_dot_freq = process_tracks(base_diname, base_star[0], i, suffix, history_finame, orbev_finame)
		t_ind=np.argmin(np.abs(time-t0s[i]))
		e0=e[t_ind]
		a0=a[t_ind]
		t0=time[t_ind]
		e_eq, a_eq = equilibrium_evolution(time, e_dot_freq, a_dot_freq, e0=e0, a0=a0, t0=t0)
		add_trajectory(axs[0], time-t0, e_eq, a_eq, c=colors[i], ls=linestyles[i], lw=linewidths[i], a_set=a_set)

	# variation in semi-major axis
	for i in range(4,7):
		time, e, a, e_dot_freq, a_dot_freq = process_tracks(base_diname, base_star[0], i, suffix, history_finame, orbev_finame)
		t_ind=np.argmin(np.abs(time-t0s[i]))
		e0=e[t_ind]
		a0=a[t_ind]
		t0=time[t_ind]
		e_eq, a_eq = equilibrium_evolution(time, e_dot_freq, a_dot_freq, e0=e0, a0=a0, t0=t0)
		add_trajectory(axs[1], time-t0, e_eq, a_eq, c=colors[i], ls=linestyles[i], lw=linewidths[i], a_set=a_set)
	add_trajectory(axs[1], time_base-t0, e_eq_base, a_eq_base, a_set=a_set)

	# variation in stellar rotation rate
	for i in range(7,10):
		time, e, a, e_dot_freq, a_dot_freq = process_tracks(base_diname, base_star[0], i, suffix, history_finame, orbev_finame)
		t_ind=np.argmin(np.abs(time-t0s[i]))
		e0=e[t_ind]
		a0=a[t_ind]
		t0=time[t_ind]
		e_eq, a_eq = equilibrium_evolution(time, e_dot_freq, a_dot_freq, e0=e0, a0=a0, t0=t0)
		add_trajectory(axs[2], time-t0, e_eq, a_eq, c=colors[i], ls=linestyles[i], lw=linewidths[i], a_set=a_set)
	add_trajectory(axs[2], time_base-t0, e_eq_base, a_eq_base, a_set=a_set)

	# variation in initial time
	for i in range(10,12):
		time, e, a, e_dot_freq, a_dot_freq = process_tracks(base_diname, base_star[0], i, suffix, history_finame, orbev_finame)
		t_ind=np.argmin(np.abs(time-t0s[i]))
		e0=e[t_ind]
		a0=a[t_ind]
		t0=time[t_ind]
		e_eq, a_eq = equilibrium_evolution(time, e_dot_freq, a_dot_freq, e0=e0, a0=a0, t0=t0)
		add_trajectory(axs[3], time-t0, e_eq, a_eq, c=colors[i], ls=linestyles[i], lw=linewidths[i], a_set=a_set)
	add_trajectory(axs[3], time_base-t0, e_eq_base, a_eq_base, a_set=a_set)

	# variation in mass
	for i in range(12,14):
		time, e, a, e_dot_freq, a_dot_freq = process_tracks(base_diname, base_star[i-11], i, suffix, history_finame, orbev_finame)
		t_ind=np.argmin(np.abs(time-t0s[i]))
		e0=e[t_ind]
		a0=a[t_ind]
		t0=time[t_ind]
		e_eq, a_eq = equilibrium_evolution(time, e_dot_freq, a_dot_freq, e0=e0, a0=a0, t0=t0)
		add_trajectory(axs[4], time-t0, e_eq, a_eq, c=colors[i], ls=linestyles[i], lw=linewidths[i], a_set=a_set)
	add_trajectory(axs[4], time_base-t0, e_eq_base, a_eq_base, a_set=a_set)

	# variation in metallicity
	for i in range(14,16):
		time, e, a, e_dot_freq, a_dot_freq = process_tracks(base_diname, base_star[i-11], i, suffix, history_finame, orbev_finame)
		t_ind=np.argmin(np.abs(time-t0s[i]))
		e0=e[t_ind]
		a0=a[t_ind]
		t0=time[t_ind]
		e_eq, a_eq = equilibrium_evolution(time, e_dot_freq, a_dot_freq, e0=e0, a0=a0, t0=t0)
		add_trajectory(axs[5], time-t0, e_eq, a_eq, c=colors[i], ls=linestyles[i], lw=linewidths[i], a_set=a_set)
	add_trajectory(axs[5], time_base-t0, e_eq_base, a_eq_base, a_set=a_set)

	plt.show()

def process_tracks(base_diname, base_star, i, suffix, history_finame, orbev_finame):
	base_finame = base_diname+"{}_orb{}_{}".format(base_star, i, suffix)
	# load base data
	fwd_fail=False
	rev_fail=False
	try:
		# load forward and reverse tracks
		time_fwd, e_fwd, a_fwd, e_dot_fwd, a_dot_fwd = process_single_track(base_finame+"_fwd/", history_finame, orbev_finame)
	except Exception as ex:
		print(ex)
		print("Failed to load forward track for {}_orb{}".format(base_star, i))
		fwd_fail=True

	try:
		time_rev, e_rev, a_rev, e_dot_rev, a_dot_rev = process_single_track(base_finame+"_rev/", history_finame, orbev_finame)
	except Exception as ex:
		print(ex)
		print("Failed to load reverse track for {}_orb{}".format(base_star, i))
		rev_fail=True

	if fwd_fail and rev_fail:
		return None, None, None, None, None
	elif fwd_fail:
		time=time_rev[::-1]
		e=e_rev[::-1]
		a=a_rev[::-1]

		e_dot=e_dot_rev[::-1]
		a_dot=a_dot_rev[::-1]
	elif rev_fail:
		time=time_fwd
		e=e_fwd
		a=a_fwd

		e_dot=e_dot_fwd
		a_dot=a_dot_fwd
	else:
		# joint forward and reverse tracks
		time=join_fwd_rev_tracks(time_rev, time_fwd)
		e=join_fwd_rev_tracks(e_rev, e_fwd)
		a=join_fwd_rev_tracks(a_rev, a_fwd)

		e_dot=join_fwd_rev_tracks(e_dot_rev, e_dot_fwd)
		a_dot=join_fwd_rev_tracks(a_dot_rev, a_dot_fwd)

	#return time, e, a, e_dot_freq, a_dot_freq
	return time, e, a, e_dot, a_dot

def load_history(finame):
	h=pd.read_csv(finame)
	return h["time"].values, h["e"].values, h["a"].values

def load_dimensions(finame):
	d=np.loadtxt(finame)
	return d

def load_response(finame):
	with h5py.File(finame, "r") as hf:
		e_dot_freq=hf["e_dot_freq"][:]
		a_dot_freq=hf["a_dot_freq"][:]

	return e_dot_freq, a_dot_freq

def dimensionalize(data, fs, Rs=None):
	data*=365
	for i in range(data.shape[0]):
		data[i]*=fs[i]

		if Rs is not None:
			data[i]*=(Rs[i]/au)

	return data

def process_single_track(fidir, history_finame, orbev_finame):
	# load data
	time, e, a = load_history(fidir+history_finame)
	e_dot=(e[1:]-e[:-1])/(time[1:]-time[:-1])
	a_dot=(a[1:]-a[:-1])/(time[1:]-time[:-1])

	good_inds=((~np.isnan(e_dot) & ~np.isnan(a_dot)) & (~np.isinf(e_dot) & ~np.isinf(a_dot))) & ((e_dot!=0) & (a_dot!=0))
	time=time[1:][good_inds]
	e=e[1:][good_inds]
	a=a[1:][good_inds]
	e_dot=e_dot[good_inds]
	a_dot=a_dot[good_inds]

	# trim
	npts=np.min([len(time), len(e_dot)])
	e=e[:npts]
	a=a[:npts]

	e_dot=e_dot[:npts]
	a_dot=a_dot[:npts]
	time=time[:npts]

	return time, e, a, e_dot, a_dot 

def join_fwd_rev_tracks(data_rev, data_fwd):
	return np.concatenate([data_rev[::-1,...], data_fwd[1:,...]])

def equilibrium_evolution(time, e_dot_freq, a_dot_freq, e0=0.5, a0=0.06, t0=0):
	# calculate total tidal evolution rate (equilibrium + dynamical)
	e_dot_total=e_dot_freq
	a_dot_total=a_dot_freq

	# initialize equilibrium rate working array
	e_dot_equilibrium=np.zeros(e_dot_freq.shape[0])
	a_dot_equilibrium=np.zeros(e_dot_freq.shape[0])

	win_len=5000
	for i in range(len(e_dot_equilibrium)):
		if i%win_len==0:
			print("{}/{}".format(i, len(e_dot_equilibrium)))

		# slice total orbital evolution rate, respecting boundaries of the array
		ind1=i-win_len
		ind2=i+win_len
		if i<win_len:
			ind1=0
		if i>(len(e_dot_equilibrium-win_len)-win_len):
			ind2=len(e_dot_equilibrium)

		avg_inds=(time>=time[i]-50e6) & (time<=time[i]+50e6)

		# define equilibrium tidal evolution rate as the 1% evolution rate over 5000 nearest samples (empirical)
		e_dot_equilibrium[i]=np.sign(e_dot_total[ind1])*np.percentile(np.abs(e_dot_total[avg_inds]), 1)
		a_dot_equilibrium[i]=np.sign(a_dot_total[ind1])*np.percentile(np.abs(a_dot_total[avg_inds]), 1)

	# initialize working arrays and set initial orbital configuration
	t_start_ind=np.argmin(np.abs(time-t0))
	# "equilibrium" values
	e_eq=np.zeros(len(time))
	e_eq[t_start_ind]=e0
	a_eq=np.zeros(len(time))
	a_eq[t_start_ind]=a0

	dt=time[1:]-time[:-1]
	dt=np.concatenate([[dt[0]], dt])

	for i in range(t_start_ind,0,-1):
		e_eq[i-1] = e_eq[i] - dt[i]*e_dot_equilibrium[i]
		a_eq[i-1] = a_eq[i] - dt[i]*a_dot_equilibrium[i]

	for i in range(t_start_ind,len(time)-1):
		e_eq[i+1] = e_eq[i] + dt[i]*e_dot_equilibrium[i]
		a_eq[i+1] = a_eq[i] + dt[i]*a_dot_equilibrium[i]

	return e_eq, a_eq

def add_trajectory(axs, time, e, a, c='k', alpha=1.0, lw=1, ls='-', a_set=0):
	t_ind=np.argmin(np.abs(time))
	a=(a-a[t_ind]) + a_set
	time=time/1e6
	axs[0].plot(time, e, c=c, alpha=alpha, lw=lw, linestyle=ls)
	axs[1].plot(time, a, c=c, alpha=alpha, lw=lw, linestyle=ls)

	axs[0].scatter(time[-1], e[-1], c=c, marker='x')
	axs[1].scatter(time[-1], a[-1], c=c, marker='x')

	axs[0].set_ylabel("e",fontsize=14)
	axs[1].set_ylabel("a [au]",fontsize=14)

	axs[0].set_xticks([-1500,-1000,-500,0,500,1000,1500])
	axs[1].set_xticks([-1500,-1000,-500,0,500,1000,1500])

if __name__=="__main__":
	main()