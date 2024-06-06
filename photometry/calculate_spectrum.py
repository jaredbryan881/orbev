import numpy as np
import pandas as pd
import h5py
import os

import pymsg as pm

import scipy.special as ss

def main():
	base_finame="../data/long_time/M1.36_Z0.02262_orb0_loworot"
	with h5py.File(base_finame+"/tidal_response_history.h5", "r") as hf:
		lag_L=hf["lag_L_ref"][:]
		xi_r=hf["xi_r_ref"][:]
		Omega_orb=hf["Omega_orb"][:]

	h=pd.read_csv("{}/orbital_history.data".format(base_finame))
	fs=np.loadtxt(base_finame+"/fs.txt")
	Teff=np.loadtxt(base_finame+"/Ts.txt")
	logg=np.loadtxt(base_finame+"/loggs.txt")

	npts=np.min([h.time.values.shape[0], len(fs)])

	n_times=npts
	n_k=lag_L.shape[2]

	theta_0=0
	phi_0=0

	# Load MSG spectroscopic grid
	MSG_DIR = os.environ['MSG_DIR']
	GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')

	specgrid_file_name = os.path.join(GRID_DIR, 'sg-demo.h5')
	specgrid = pm.SpecGrid(specgrid_file_name)

	# define wavelength bands
	n_lam=7
	lams=np.linspace(6000., 9000., n_lam) # [Angstrom]
	n_bands=n_lam-1

	# initialize the amplitude spectra
	dJ = np.zeros((n_times,n_k,n_bands), dtype=complex)
	# iterate over times
	for i_t in range(n_times):
		# calculate bandpass correction for current temperature
		atm_params = {'Teff': Teff[i_t], 'log(g)': logg[i_t]}
		omega = np.arange(n_k)*Omega_orb[i_t]/fs[i_t]

		# calculate amplitude spectrum for each wavelength band and time
		try:
			dJ[i_t]=calculate_spectrum(xi_r[i_t], lag_L[i_t], omega, specgrid, lams, atm_params, theta_0, phi_0)
		except Exception as ex:
			print(ex)
			dJ[i_t]=np.nan

	time=h.time.values[:npts]
	dJ=dJ[:npts,:,:]

	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	for i_b in range(n_bands):
		dJ_net=dJ[:,:,i_b].sum(axis=1)
		plt.plot(time, np.log10(np.abs(dJ_net)), c=cm.inferno(i_b/n_bands))
	plt.show()

	for i_b in range(n_bands):
		plt.imshow(np.log10(np.abs(dJ[:,:,i_b].T)), origin='lower', cmap='inferno', clim=[-9,-3])
		plt.gca().set_aspect("auto")
		plt.show()

def calculate_spectrum(xi_r, lag_L, omega, specgrid, lams, x, theta_0, phi_0):
	"""
	Evaluate Fourier amplitudes.
	"""
	# Initialize frequency / amplitude arrays

	n_k=xi_r.shape[1]
	n_lam=len(lams)

	l = 2

	ks = np.arange(n_k)
	A = np.zeros((n_k,n_lam-1), dtype=complex)

	# Iterate through rows
	for (i_m,m) in enumerate([0,2]):
		# Townsend (2003), eqns. 4-6
		Del_R = xi_r[i_m]
		Del_L = lag_L[i_m]
		Del_T = 0.25*(Del_L - 2*Del_R)
		Del_g = -(2 + omega**2)*Del_R

		# Townsend (2003), eqn. 15
		I_0 = specgrid.D_moment(x, 0, lams)
		I_l = specgrid.D_moment(x, l, lams)

		dI_l_dlnTeff = specgrid.D_moment(x, l, lams, deriv={'Teff': True})*x['Teff']
		dI_l_dlng    = specgrid.D_moment(x, l, lams, deriv={'log(g)': True})/np.log(10)

		# Townsend (2003), eqns. 12-14
		Y_lm = ss.sph_harm(m, l, phi_0, theta_0)

		R_lm = (2+l)*(1-l)*I_l/I_0*Y_lm
		T_lm = dI_l_dlnTeff/I_0*Y_lm
		G_lm = dI_l_dlng/I_0*Y_lm

		# Willems+(2010), eqn. 53
		kappas=get_kappa(ks, m)

		for i_b in range(n_lam-1):
			# Townsend (2003), eqn. 11
			A[:,i_b] += 2*kappas*(Del_R*R_lm[i_b] + Del_T*T_lm[i_b] + Del_g*G_lm[i_b])

	# Return data
	return A

def get_kappa(ks, m):
	kappas=[]
	for k in ks:
		if k == 0:
			if m == 0:
				kappa = 0.5
			elif m >= 1:
				kappa = 1.
			else:
				kappa = 0.
		else:
			kappa = 1.
		kappas.append(kappa)

	return np.array(kappas)

if __name__=="__main__":
	main()