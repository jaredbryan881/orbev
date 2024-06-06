import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import mesa_reader as mr
import os

import warnings
warnings.filterwarnings("ignore")

def main():
	base_fidir="/home/jaredbryan/MIT/astero/gyre_HATP2/orbev/data/long_time/"
	base_finame="M1.36_Z0.02262_orb0_highorot"
	#foname="edot.png"
	plot_orbev_rates(base_fidir+base_finame, foname=None)

Rsun = 695700000 # [m]
au = 1.496e+11   # [m]

def plot_orbev_rates(fidir, foname=None):
	orbev_finame="tidal_response_history.h5"
	history_finame="orbital_history.data"

	# load forward data
	with h5py.File(fidir+"_fwd/"+orbev_finame, "r") as hf:
		a_dot_freq = hf["a_dot_freq"][:]
		e_dot_freq = hf["e_dot_freq"][:]
		o_dot_freq = hf["o_dot_freq"][:]
		J_dot_freq = hf["J_dot_freq"][:]
	data=pd.read_csv(fidir+"_fwd/"+history_finame)
	time=data["time"]/1e6

	# trim
	npts=np.min([len(time), a_dot_freq.shape[0]])
	a_dot_freq=a_dot_freq[:npts,:,:]
	o_dot_freq=o_dot_freq[:npts,:,:]
	e_dot_freq=e_dot_freq[:npts,:,:]
	J_dot_freq=J_dot_freq[:npts,:,:]
	time=time[:npts]

	# dimensionalize
	fs=np.loadtxt(fidir+"_fwd/fs.txt")
	Rs=np.loadtxt(fidir+"_fwd/Rs.txt")
	orbev_freq = e_dot_freq*365
	for i in range(orbev_freq.shape[0]):
		orbev_freq[i,:,:]*=fs[i]
		#orbev_freq[i,:,:]*=(Rs[i]/au)

	# load reverse data
	with h5py.File(fidir+"_rev/"+orbev_finame, "r") as hf:
		a_dot_freq_rev = hf["a_dot_freq"][:]
		e_dot_freq_rev = hf["e_dot_freq"][:]
		o_dot_freq_rev = hf["o_dot_freq"][:]
		J_dot_freq_rev = hf["J_dot_freq"][:]
	data_rev=pd.read_csv(fidir+"_rev/"+history_finame)
	time_rev=data_rev["time"]/1e6

	fs=np.loadtxt(fidir+"_rev/"+"fs.txt")
	Rs=np.loadtxt(fidir+"_rev/Rs.txt")

	# trim
	npts=np.min([len(time), a_dot_freq_rev.shape[0], len(fs)])
	a_dot_freq_rev=a_dot_freq_rev[:npts,:,:]
	o_dot_freq_rev=o_dot_freq_rev[:npts,:,:]
	e_dot_freq_rev=e_dot_freq_rev[:npts,:,:]
	J_dot_freq_rev=J_dot_freq_rev[:npts,:,:]
	time_rev=time_rev[:npts]

	# dimensionalize
	orbev_freq_rev = e_dot_freq_rev*365
	for i in range(orbev_freq_rev.shape[0]):
		orbev_freq_rev[i,:,:]*=fs[i]
		#orbev_freq_rev[i,:,:]*=(Rs[i]/au)

	time=np.concatenate((time_rev.values[::-1], time.values))
	orbev_freq = np.concatenate((orbev_freq_rev[::-1,:,:], orbev_freq), axis=0)

	grid_time, grid_k = np.meshgrid(time, np.arange(orbev_freq.shape[2]))

	clim=[-20,-6]

	fig1,axs1=plt.subplots(2,2,sharex=True,sharey='row',figsize=(10,5))

	# m=0
	orbev_dot_mm2=orbev_freq[:,0,:]
	orbev_dot_mm2_net=np.sum(orbev_dot_mm2, axis=1)
	axs1[0,0].scatter(time[orbev_dot_mm2_net>0], np.log10(np.abs(orbev_dot_mm2_net[orbev_dot_mm2_net>0])), color='crimson', s=1)
	axs1[0,0].scatter(time[orbev_dot_mm2_net<0], np.log10(np.abs(orbev_dot_mm2_net[orbev_dot_mm2_net<0])), color='steelblue', s=1)
	orbev_dot_mm2_pos=np.ma.masked_where(orbev_dot_mm2<0, orbev_dot_mm2)
	orbev_dot_mm2_neg=np.ma.masked_where(orbev_dot_mm2>0, orbev_dot_mm2)
	axs1[1,0].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_mm2_pos)).T, cmap='Reds', vmin=clim[0], vmax=clim[1], shading='auto')
	axs1[1,0].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_mm2_neg)).T, cmap='Blues', vmin=clim[0], vmax=clim[1], shading='auto')

	# m=2
	orbev_dot_m0=orbev_freq[:,1,:]
	orbev_dot_m0_net=np.sum(orbev_dot_m0, axis=1)
	axs1[0,1].scatter(time[orbev_dot_m0_net>0], np.log10(np.abs(orbev_dot_m0_net[orbev_dot_m0_net>0])), color='crimson', s=1)
	axs1[0,1].scatter(time[orbev_dot_m0_net<0], np.log10(np.abs(orbev_dot_m0_net[orbev_dot_m0_net<0])), color='steelblue', s=1)
	orbev_dot_m0_pos=np.ma.masked_where(orbev_dot_m0<0, orbev_dot_m0)
	orbev_dot_m0_neg=np.ma.masked_where(orbev_dot_m0>0, orbev_dot_m0)
	axs1[1,1].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_m0_pos)).T, cmap='Reds', vmin=clim[0], vmax=clim[1], shading='auto')
	axs1[1,1].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_m0_neg)).T, cmap='Blues', vmin=clim[0], vmax=clim[1], shading='auto')

	axs1[0,0].set_ylim(-14,-3)

	axs1[0,0].set_title("m=0",fontsize=14)
	axs1[0,1].set_title("m=2",fontsize=14)

	axs1[1,0].set_xlabel("Time [Myr]",fontsize=14)
	axs1[1,1].set_xlabel("Time [Myr]",fontsize=14)

	axs1[0,0].set_ylabel(r"$\dot{e}$ [1/yr]", fontsize=14)
	axs1[1,0].set_ylabel(r"$k=\sigma_{m,k}/\Omega_{orb}$", fontsize=14)

	plt.tight_layout()
	if foname is not None:
		plt.savefig(foname)

	plt.show()

if __name__=="__main__":
	main()