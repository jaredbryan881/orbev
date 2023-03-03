import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

def main():
	for i in range(0,12):
		fidir="/home/jared/MIT/astero/gyre_HATP2/orbev/data/live_fwd_ps_orot/M1.36_Z0.02262_orb{}_live_fwd/".format(i)
		plot_orbev_rates(fidir)

def plot_orbev_rates(fidir):
	orbev_finame="tidal_response_history.h5"
	history_finame="orbital_history.data"

	with h5py.File(fidir+orbev_finame, "r") as hf:
		a_dot_freq = hf["a_dot_freq"][:]
		e_dot_freq = hf["e_dot_freq"][:]
		J_dot_freq = hf["J_dot_freq"][:]

	data=pd.read_csv(fidir+history_finame)
	time=data["time"]/1e6

	clim=[-20,-10]

	npts=np.min([len(time), a_dot_freq.shape[0]])

	a_dot_freq=a_dot_freq[:npts,:,:]
	e_dot_freq=e_dot_freq[:npts,:,:]
	J_dot_freq=J_dot_freq[:npts,:,:]
	time=time[:npts]

	orbev_freq = e_dot_freq

	grid_time, grid_k = np.meshgrid(time, np.arange(orbev_freq.shape[2]))

	fig1,axs1=plt.subplots(2,3,sharex=True,sharey='row', figsize=(10,5))

	# m=-2
	orbev_dot_mm2=orbev_freq[:,0,:]
	orbev_dot_mm2_net=np.sum(orbev_dot_mm2, axis=1)
	axs1[0,0].scatter(time[orbev_dot_mm2_net>0], np.log10(np.abs(orbev_dot_mm2_net[orbev_dot_mm2_net>0])), color='crimson', s=1)
	axs1[0,0].scatter(time[orbev_dot_mm2_net<0], np.log10(np.abs(orbev_dot_mm2_net[orbev_dot_mm2_net<0])), color='steelblue', s=1)
	orbev_dot_mm2_pos=np.ma.masked_where(orbev_dot_mm2<0, orbev_dot_mm2)
	orbev_dot_mm2_neg=np.ma.masked_where(orbev_dot_mm2>0, orbev_dot_mm2)
	axs1[1,0].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_mm2_pos)).T, cmap='Reds', vmin=clim[0], vmax=clim[1], shading='auto')
	axs1[1,0].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_mm2_neg)).T, cmap='Blues', vmin=clim[0], vmax=clim[1], shading='auto')

	# m=0
	orbev_dot_m0=orbev_freq[:,2,:]
	orbev_dot_m0_net=np.sum(orbev_dot_m0, axis=1)
	axs1[0,1].scatter(time[orbev_dot_m0_net>0], np.log10(np.abs(orbev_dot_m0_net[orbev_dot_m0_net>0])), color='crimson', s=1)
	axs1[0,1].scatter(time[orbev_dot_m0_net<0], np.log10(np.abs(orbev_dot_m0_net[orbev_dot_m0_net<0])), color='steelblue', s=1)
	orbev_dot_m0_pos=np.ma.masked_where(orbev_dot_m0<0, orbev_dot_m0)
	orbev_dot_m0_neg=np.ma.masked_where(orbev_dot_m0>0, orbev_dot_m0)
	axs1[1,1].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_m0_pos)).T, cmap='Reds', vmin=clim[0], vmax=clim[1], shading='auto')
	axs1[1,1].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_m0_neg)).T, cmap='Blues', vmin=clim[0], vmax=clim[1], shading='auto')

	# m=2
	orbev_dot_mp2=orbev_freq[:,4,:]
	orbev_dot_mp2_net=np.sum(orbev_dot_mp2, axis=1)
	axs1[0,2].scatter(time[orbev_dot_mp2_net>0], np.log10(np.abs(orbev_dot_mp2_net[orbev_dot_mp2_net>0])), color='crimson', s=1)
	axs1[0,2].scatter(time[orbev_dot_mp2_net<0], np.log10(np.abs(orbev_dot_mp2_net[orbev_dot_mp2_net<0])), color='steelblue', s=1)
	orbev_dot_mp2_pos=np.ma.masked_where(orbev_dot_mp2<0, orbev_dot_mp2)
	orbev_dot_mp2_neg=np.ma.masked_where(orbev_dot_mp2>0, orbev_dot_mp2)
	axs1[1,2].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_mp2_pos)).T, cmap='Reds', vmin=clim[0], vmax=clim[1], shading='auto')
	axs1[1,2].pcolormesh(grid_time, grid_k, np.log10(np.abs(orbev_dot_mp2_neg)).T, cmap='Blues', vmin=clim[0], vmax=clim[1], shading='auto')

	axs1[1,0].set_xlabel("Time [Myr]", fontsize=12)
	axs1[1,1].set_xlabel("Time [Myr]", fontsize=12)
	axs1[1,2].set_xlabel("Time [Myr]", fontsize=12)

	axs1[0,0].set_ylabel(r"$\log_{10}(\dot{e})$", fontsize=12)
	axs1[1,0].set_ylabel("k", fontsize=12)

	axs1[0,0].set_title("m=-2", fontsize=12)
	axs1[0,1].set_title("m=0", fontsize=12)
	axs1[0,2].set_title("m=2", fontsize=12)

	#axs1[0,0].set_xlim(46.8, 47.0)
	#axs1[0,0].set_ylim(-16,-11.5)

	plt.tight_layout()

	plt.show()

if __name__=="__main__":
	main()