import h5py
import numpy as np
import matplotlib.pyplot as plt

def main():
	finame="../data/M1.36_Z0.027_mesh0.3_gdt1e3/tidal_response_history.h5"
	with h5py.File(finame, "r") as hf:
		# harmonic degree
		l=hf["l"][:]
		# azimuthal order
		m=hf["m"][:]
		# Fourier harmonic
		k=hf["k"][:]

		# ratio of secondary mass to primary mass
		q=hf["q"][:]
		# ratio of primary radius to orbital semi-major axis
		R_a=hf["R_a"][:]
		# orbital frequency
		Omega_orb=hf["Omega_orb"][:]
		# orbital eccentricity
		e=hf["e"][:]

		# secular evolution coefficients (G factors) (see Willems et al., 2010)
		Gbar_1=hf["Gbar_1"][:]
		Gbar_2=hf["Gbar_2"][:]
		Gbar_3=hf["Gbar_3"][:]
		Gbar_4=hf["Gbar_4"][:]

		# tidal expansion coefficient (eqn. A1 of Sun et al., 2023)
		cbar=hf["cbar"][:]
		# Eulerian total potential perturbation at reference location
		eul_Psi_ref=hf["eul_Psi_ref"][:]
		# Eulerian potential perturbation at reference location
		eul_Phi_ref=hf["eul_Phi_ref"][:]
		# tidal potential at reference location
		Phi_T_ref=hf["Phi_T_ref"][:]

		# radial displacement perturbation at reference location
		xi_r_ref=hf["xi_r_ref"][:]
		# Lagrangian radiative luminosity perturbation at reference location
		lag_L_ref=hf["lag_L_ref"][:]

	# orbital parameters
	eps_tide=(R_a**3)*q

	# calculate eulerian perturbation to the gravitational potential
	eul_phi=eul_Psi_ref
	eul_phi.real-=Phi_T_ref

	# calculate tidal response amplitude and phase
	F=np.zeros(eul_Psi_ref.shape, dtype=complex)
	for t in range(F.shape[0]):
		for m in range(F.shape[1]):
			for k in range(F.shape[2]):
				F[t,m,k]=-(1/2)*(np.sqrt(4*np.pi)*eul_phi[t,m,k]/(eps_tide[t]*cbar[t,m,k]) + 1)
	gamma = np.arctan2(np.imag(F), np.real(F))

	plt.imshow(np.log10(np.abs(np.abs(F)*np.sin(gamma)*Gbar_2))[:,4,:].T, origin='lower', cmap='Blues', clim=[-20,-10])
	plt.gca().set_aspect("auto")
	plt.show()

if __name__=="__main__":
	main()