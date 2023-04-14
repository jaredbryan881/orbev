import numpy as np
import pandas as pd
import h5py

import scipy.special as ss

from bandpass_correction import bandpass_correction

def main():
	base_finame="../data/figure/gyre_maxerror/M1.36_Z0.02262_orb0_highorot"
	with h5py.File(base_finame+"/tidal_response_history.h5", "r") as hf:
		lag_L=hf["lag_L_ref"][:]
		xi_r=hf["xi_r_ref"][:]
	h=pd.read_csv("{}/orbital_history.data".format(base_finame))
	Ts=np.loadtxt(base_finame+"/Ts.txt")

	n_times=lag_L.shape[0]
	n_k=lag_L.shape[2]

	# initialize the amplitude spectra
	dJ = np.zeros((n_times,n_k), dtype=complex)
	# iterate over times
	for t in range(n_times):
		try:
			# calculate bandpass correction for current temperature
			bpc=bandpass_correction([0.5e-6,5e-6], Ts[t])

			# calculate amplitude spectrum for each time
			dJ[t]=calculate_spectrum(xi_r[t], lag_L[t], 0, 0)*bpc
		except:
			dJ[t]=np.nan

	print(np.max(np.log10(np.abs(dJ[~np.isnan(dJ)]))))

def calculate_spectrum(xi_r, lag_L, theta_0, phi_0):
	"""
	Calculate the amplitude spectrum of relative flux variations due to tidal forcing.

	This result is from Townsend (2003) A semi-analytical formula for the light variations due to low-frequency g modes in rotating stars
	This implementation is from GYRE https://github.com/rhdtownsend/gyre/blob/c1b346dabfd6ec3b30b04f06ddbc90d164fabdc8/src/tide/plot_lightcurve.py

	Arguments
	---------
	:param xi_r: np.array
		Radial displacement at the stellar surface. Shape is (n_m, n_k) and elements are complex.
	:param lag_L: np.array
		Lagrangian flux perturbation at the stellar surface. Shape is (n_m, n_k) and elements are complex.
	:param theta_0: float
		Observer coordinate
	:param phi_0: float
		Observer coordinate

	Returns
	-------
	:return A: np.array
		Amplitude spectrum of relative flux variations. Shape is (n_k) and elements are complex.
	"""
	# Initialize the amplitude spectrum for current time
	n_k=lag_L.shape[1]
	A = np.zeros(n_k, dtype=complex)

	# Townsend (2003), eqns. (15) & (16) assuming I_x = const. (no limb darkening)
	I_0 = 1/2
	I_l = 1/8 # l=2

	# Loop over m and k
	l=2
	for (i_m, m) in enumerate([0,2]):
		# Townsend (2003), eqns. (12) & (13) assuming
		# dlnI_x/dlnTeff = 4 (black body)
		Yml = ss.sph_harm(m, l, phi_0, theta_0)
		Rml = (2 + l)*(1 - l)*(I_l/I_0)*Yml
		Tml = 4*(I_l/I_0)*Yml

		# Townsend (2003), eqns. 17 & 18
		Del_R = np.sqrt(4.*np.pi)*xi_r[i_m]
		Del_L = np.sqrt(4.*np.pi)*lag_L[i_m]
		Del_T = (1/4)*(Del_L - 2*Del_R)

		# Add the Fourier contribution
		A += Del_R*Rml + Del_T*Tml

	return A

if __name__=="__main__":
	main()