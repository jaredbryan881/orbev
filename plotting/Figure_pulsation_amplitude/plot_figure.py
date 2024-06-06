import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from cycler import cycler
import mesa_reader as mr
import os

import pickle as pkl

import pymsg as pm

import sys
sys.path.append("../../photometry/")
from calculate_spectrum import calculate_spectrum

import warnings
warnings.filterwarnings("ignore")

def main():
	#base_fidir="/home/jared/MIT/astero/gyre_HATP2/orbev/data/figure/live_planet/"
	base_fidir="/home/jared/MIT/astero/gyre_HATP2/orbev/data/figure/gyre_maxerror/"
	#calc_response_evolution(base_fidir,"M1.36_Z0.02262_orb0_highorot")

	grid_time, grid_k, lams, dJ_full = load_response_evolution()

	plot_response_evolution(grid_time, grid_k, lams, dJ_full, clim=[-5,-3])

	for t_ind in range(1000,60000,50):
		print(grid_time[0,t_ind]/1e9)
		if grid_time[0,t_ind]/1e9<2.586: continue
		plot_response_spectrum(dJ_full, t_ind=t_ind)
		plot_lightcurve(dJ_full, t_ind=t_ind)

def save_response_evolution(grid_time, grid_k, lams, dJ_full):
	with open("dJ_lams.pkl", "wb") as f:
		pkl.dump([grid_time, grid_k, lams, dJ_full], f)

def load_response_evolution():
	with open("dJ_lams.pkl", "rb") as f:
		data=pkl.load(f)
	grid_time, grid_k, lams, dJ_full = data[0],data[1],data[2],data[3]
	return grid_time, grid_k, lams, dJ_full

def calc_response_evolution(base_fidir, core_fidir):
	orbev_finame="tidal_response_history.h5"
	history_finame="orbital_history.data"

	theta_0=np.pi/2
	phi_0=np.pi/4

	# Load MSG spectroscopic grid
	MSG_DIR = os.environ['MSG_DIR']
	GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')
	specgrid_file_name="sg-CAP18-coarse.h5"
	specgrid = pm.SpecGrid(specgrid_file_name)

	# define wavelength bands
	n_lam=7
	lams=np.linspace(6e3, 5e4, n_lam) # [Angstrom]
	n_bands=n_lam-1

	# load forward data
	fidir=base_fidir+core_fidir+"_fwd"
	with h5py.File(fidir+"/"+orbev_finame, "r") as hf:
		lag_L = hf["lag_L_ref"][:]
		xi_r = hf["xi_r_ref"][:]
		Omega_orb=hf["Omega_orb"][:]

	data=pd.read_csv(fidir+"/"+history_finame)
	fs=np.loadtxt(fidir+"/fs.txt")
	Teff=np.loadtxt(fidir+"/Ts.txt")
	logg=np.loadtxt(fidir+"/loggs.txt")

	npts=np.min([data.time.values.shape[0], len(fs)])

	time=data["time"]/1e6
	n_times=npts
	n_k=lag_L.shape[2]

	# initialize the amplitude spectra
	dJ = np.zeros((n_times,n_k,n_bands), dtype=complex)
	# iterate over times
	for i_t in range(n_times):
		if i_t%1000==0:
			print(i_t)
		# calculate bandpass correction for current temperature
		atm_params = {'Teff': Teff[i_t], 'log(g)': logg[i_t], '[Fe/H]': 0.02262}
		omega = np.arange(n_k)*Omega_orb[i_t]

		# calculate amplitude spectrum for each wavelength band and time
		try:
			dJ[i_t]=calculate_spectrum(xi_r[i_t], lag_L[i_t], omega, specgrid, lams, atm_params, theta_0, phi_0)
		except Exception as ex:
			print(ex)
			dJ[i_t]=np.nan

	# trim
	time=data.time.values[:npts]
	dJ=dJ[:npts,:,:]

	# load reverse data
	fidir_rev=base_fidir+core_fidir+"_rev"
	with h5py.File(fidir_rev+"/"+orbev_finame, "r") as hf:
		lag_L_rev = hf["lag_L_ref"][:]
		xi_r_rev = hf["xi_r_ref"][:]
		Omega_orb = hf["Omega_orb"][:]

	data_rev=pd.read_csv(fidir_rev+"/"+history_finame)
	fs_rev=np.loadtxt(fidir_rev+"/fs.txt")
	Teff_rev=np.loadtxt(fidir_rev+"/Ts.txt")
	logg_rev=np.loadtxt(fidir_rev+"/loggs.txt")

	npts_rev=np.min([data_rev.time.values.shape[0], len(fs_rev)])

	time_rev=data_rev["time"]/1e6
	n_times=npts_rev
	n_k=lag_L.shape[2]

	# initialize the amplitude spectra
	dJ_rev = np.zeros((n_times,n_k,n_bands), dtype=complex)
	# iterate over times
	for i_t in range(n_times):
		# calculate bandpass correction for current temperature
		atm_params = {'Teff': Teff_rev[i_t], 'log(g)': logg_rev[i_t], '[Fe/H]': 0.02262}
		omega = np.arange(n_k)*Omega_orb[i_t]

		# calculate amplitude spectrum for each wavelength band and time
		try:
			dJ_rev[i_t]=calculate_spectrum(xi_r_rev[i_t], lag_L_rev[i_t], omega, specgrid, lams, atm_params, theta_0, phi_0)
		except Exception as ex:
			dJ_rev[i_t]=np.nan

	# trim
	time_rev=data_rev.time.values[:npts_rev]
	dJ_rev=dJ_rev[:npts_rev,:,:]

	grid_time, grid_k = np.meshgrid(time, np.arange(dJ.shape[1]))

	time_full = np.concatenate((time_rev[::-1], time))
	dJ_full = np.concatenate((dJ_rev[::-1,:], dJ), axis=0)

	print(np.max(np.log10(np.abs(dJ_full))))

	grid_time, grid_k = np.meshgrid(time_full, np.arange(dJ_full.shape[1]))

	save_response_evolution(grid_time, grid_k, lams, dJ_full)

def plot_response_spectrum(A, t_ind=1000):
	ks = np.arange(A.shape[1])
	fig,ax=plt.subplots(1,1)
	for i_b in range(A.shape[2]):
		dk=i_b/(A.shape[2])
		ax.vlines(x=ks+dk, ymin=-14, ymax=np.log10(np.abs(A[int(t_ind),:,i_b])), colors=cm.inferno(i_b/(1.5*A.shape[2])))
	ax.set_ylabel(r'$log_{10}(|dF/F|)$', fontsize=14)
	ax.set_xlabel(r"$k=\sigma_{m,k}/\Omega_{orb}$", fontsize=14)
	ax.set_ylim(-5,-4)
	ax.set_xlim(0,50)
	plt.tight_layout()
	plt.show()

def gt(A_cur, A_med, factor=10):
	return np.any(A_cur>factor*A_med)

def plot_response_evolution(grid_time, grid_k, lams, A, l_ind=None, clim=[-9,-3]):
	if l_ind is None:
		A=A.sum(axis=2)
	else:
		A=A[:,:,l_ind]

	times=grid_time[0,:]
	indicator=np.sum(np.abs(A)>100*np.median(np.abs(A),axis=0), axis=1)

	plt.plot(times, indicator)
	plt.show()

	plt.plot(np.max(np.log10(np.abs(A)), axis=0))
	plt.plot(np.min(np.log10(np.abs(A)), axis=0))
	plt.plot(np.median(np.log10(np.abs(A)), axis=0))
	plt.show()

	fig,axs=plt.subplots(1,1,figsize=(10,5))
	axs.pcolormesh(grid_time/1e6, grid_k, (np.log10(np.abs(A))).T, cmap='inferno', vmin=clim[0], vmax=clim[1], shading='auto')
	axs.set_ylabel(r"$k=\sigma_{m,k}/\Omega_{orb}$", fontsize=14)
	axs.set_xlabel("Time [Myr]", fontsize=14)
	plt.tight_layout()
	plt.show()

def plot_lightcurve(A, t_ind=1000):
	# Evaluate the light curves
	n_phase = 5001
	n_freq  = A.shape[0]
	n_k     = A.shape[1]
	n_bands = A.shape[2]

	ks=np.arange(n_k)

	phase = np.linspace(0, 4*np.pi, n_phase)
	dF = np.zeros((n_phase,n_bands))

	for i in range(n_phase):
		expp = np.exp(1j*ks*phase[i])
		dF[i,:] += np.real(np.matmul(expp, A[t_ind,:,:]))

	fig,ax=plt.subplots(1,1)
	for j in range(n_bands):
		ax.plot(phase/(2*np.pi), dF[:,j], c=cm.inferno(j/(1.5*n_bands)))

	plt.xlabel('phase/2pi',fontsize=20)
	plt.ylabel('dF/F',fontsize=20)
	plt.ylim(-0.00015, 0.00015)
	plt.xlim(0,2)
	plt.show()

if __name__=="__main__":
	main()