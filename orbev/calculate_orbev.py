import os
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

def main():
	# load command line arguments
	orbev_fodir=sys.argv[1]
	pind=int(sys.argv[2])

	try:
		with h5py.File("{}/tide_orbit.h5".format(orbev_fodir), "r") as hf:
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
	except OSError as ex:
		# provide exit code to stop simulation
		print("Failed to load tide_orbit.h5")
		sys.exit(1)

	# orbital parameters
	Omega_orb=Omega_orb[0]
	R_a=R_a[0]
	q=q[0]
	e=e[0]
	eps_tide=(R_a**3)*q

	# calculate eulerian perturbation to the gravitational potential
	eul_Psi_ref_temp=np.zeros(len(eul_Psi_ref), dtype=np.complex128)
	for i in range(len(eul_Psi_ref)):
		eul_Psi_ref_temp[i] = eul_Psi_ref[i][0] + 1j*eul_Psi_ref[i][1]
	eul_phi=eul_Psi_ref_temp
	eul_phi.real-=Phi_T_ref

	lun=np.unique(l)
	mun=np.unique(m)
	kun=np.unique(k)

	n_l=len(lun)
	n_m=len(mun)
	n_k=len(kun)

	# orbev rates
	o_dots=np.zeros((n_l,n_m,n_k))
	a_dots=np.zeros((n_l,n_m,n_k))
	e_dots=np.zeros((n_l,n_m,n_k))
	J_dots=np.zeros((n_l,n_m,n_k))

	# fundamental parameters used to calculate orbev rates
	Gbar_1_temp=np.zeros((n_l,n_m,n_k))
	Gbar_2_temp=np.zeros((n_l,n_m,n_k))
	Gbar_3_temp=np.zeros((n_l,n_m,n_k))
	Gbar_4_temp=np.zeros((n_l,n_m,n_k))
	cbar_temp=np.zeros((n_l,n_m,n_k))
	eul_Phi_ref_temp=np.zeros((n_l,n_m,n_k), dtype=complex)
	eul_Psi_ref_temp=np.zeros((n_l,n_m,n_k), dtype=complex)
	Phi_T_ref_temp=np.zeros((n_l,n_m,n_k))

	# surface response
	xi_r_ref_temp=np.zeros((n_l,n_m,n_k), dtype=complex)
	lag_L_ref_temp=np.zeros((n_l,n_m,n_k), dtype=complex)

	for ind in range(len(l)):
		# locate the current mode
		i_l = np.where(lun==l[ind])[0]
		i_m = np.where(mun==m[ind])[0]
		i_k = np.where(kun==k[ind])[0]

		c=cbar[ind]
		kappa=get_kappa(m[ind], k[ind])
		if (c==0) or (kappa==0):
			continue

		# calculate tidal response amplitude and phase
		F = -(1/2)*(np.sqrt(4*np.pi)*eul_phi[ind]/(eps_tide*c) + 1)
		gamma = np.arctan2(np.imag(F), np.real(F))

		# calculate orbital evolution rates
		o_dots[i_l,i_m,i_k] = 4*Omega_orb*q*(R_a)**(l[ind]+3)*kappa*np.abs(F)*np.cos(gamma)*Gbar_1[ind]
		a_dots[i_l,i_m,i_k] = 4*Omega_orb*(q/R_a)*(R_a)**(l[ind]+3)*kappa*np.abs(F)*np.sin(gamma)*Gbar_2[ind]
		e_dots[i_l,i_m,i_k] = 4*Omega_orb*q*(R_a)**(l[ind]+3)*kappa*np.abs(F)*np.sin(gamma)*Gbar_3[ind]
		J_dots[i_l,i_m,i_k] = 4*Omega_orb*q**2/np.sqrt(R_a*(1+q))*(R_a)**(l[ind]+3)*kappa*np.abs(F)*np.sin(gamma)*Gbar_4[ind]

		Gbar_1_temp[i_l,i_m,i_k] = Gbar_1[ind]
		Gbar_2_temp[i_l,i_m,i_k] = Gbar_2[ind]
		Gbar_3_temp[i_l,i_m,i_k] = Gbar_3[ind]
		Gbar_4_temp[i_l,i_m,i_k] = Gbar_4[ind]
		cbar_temp[i_l,i_m,i_k] = cbar[ind]
		
		eul_Phi_ref_temp[i_l,i_m,i_k] = eul_Phi_ref[ind][0]+eul_Phi_ref[ind][1]*1j
		eul_Psi_ref_temp[i_l,i_m,i_k] = eul_Psi_ref[ind][0]+eul_Psi_ref[ind][1]*1j
		Phi_T_ref_temp[i_l,i_m,i_k] = Phi_T_ref[ind]

		xi_r_ref_temp[i_l,i_m,i_k] = xi_r_ref[ind][0]+xi_r_ref[ind][1]*1j
		lag_L_ref_temp[i_l,i_m,i_k] = lag_L_ref[ind][0]+lag_L_ref[ind][1]*1j

	Gbar_1=Gbar_1_temp
	Gbar_2=Gbar_2_temp
	Gbar_3=Gbar_3_temp
	Gbar_4=Gbar_4_temp
	cbar=cbar_temp
	eul_Phi_ref=eul_Phi_ref_temp
	eul_Psi_ref=eul_Psi_ref_temp
	Phi_T_ref=Phi_T_ref_temp
	xi_r_ref=xi_r_ref_temp
	lag_L_ref=lag_L_ref_temp

	if os.path.exists("{}/tidal_response_history.h5".format(orbev_fodir)):
		# Resize datasets and append current values
		with h5py.File("{}/tidal_response_history.h5".format(orbev_fodir), "a") as hf:
			# Append current orbital parameters
			# ratio of secondary mass to primary mass
			hf["q"].resize((hf["q"].shape[0]+1), axis=0)
			hf["q"][-1]=np.array([q])
			# ratio of primary radius to semi-major axis
			hf["R_a"].resize((hf["R_a"].shape[0]+1), axis=0)
			hf["R_a"][-1]=np.array([R_a])
			# orbital frequency
			hf["Omega_orb"].resize((hf["Omega_orb"].shape[0]+1), axis=0)
			hf["Omega_orb"][-1]=np.array([Omega_orb])
			# eccentricity
			hf["e"].resize((hf["e"].shape[0]+1), axis=0)
			hf["e"][-1]=np.array([e])

			# Append rate-changes in orbital parameters
			# argument of periastron
			hf["o_dot_freq"].resize((hf["o_dot_freq"].shape[0]+1), axis=0)
			hf["o_dot_freq"][-1]=o_dots
			hf["o_dot"].resize((hf["o_dot"].shape[0]+1), axis=0)
			hf["o_dot"][-1]=np.array([np.sum(o_dots)])
			# semi-major axis
			hf["a_dot_freq"].resize((hf["a_dot_freq"].shape[0]+1), axis=0)
			hf["a_dot_freq"][-1]=a_dots
			hf["a_dot"].resize((hf["a_dot"].shape[0]+1), axis=0)
			hf["a_dot"][-1]=np.array([np.sum(a_dots)])
			# eccentricity
			hf["e_dot_freq"].resize((hf["e_dot_freq"].shape[0]+1), axis=0)
			hf["e_dot_freq"][-1]=e_dots
			hf["e_dot"].resize((hf["e_dot"].shape[0]+1), axis=0)
			hf["e_dot"][-1]=np.array([np.sum(e_dots)])
			# rotational angular momentum
			hf["J_dot_freq"].resize((hf["J_dot_freq"].shape[0]+1), axis=0)
			hf["J_dot_freq"][-1]=J_dots
			hf["J_dot"].resize((hf["J_dot"].shape[0]+1), axis=0)
			hf["J_dot"][-1]=np.array([np.sum(J_dots)])

			# Append the fundamental parameters used to calculate the rate-changes in orbital parameters
			hf["Gbar_1"].resize((hf["Gbar_1"].shape[0]+1), axis=0)
			hf["Gbar_1"][-1]=Gbar_1
			hf["Gbar_2"].resize((hf["Gbar_2"].shape[0]+1), axis=0)
			hf["Gbar_2"][-1]=Gbar_2
			hf["Gbar_3"].resize((hf["Gbar_3"].shape[0]+1), axis=0)
			hf["Gbar_3"][-1]=Gbar_3
			hf["Gbar_4"].resize((hf["Gbar_4"].shape[0]+1), axis=0)
			hf["Gbar_4"][-1]=Gbar_4
			hf["cbar"].resize((hf["cbar"].shape[0]+1), axis=0)
			hf["cbar"][-1]=cbar
			hf["eul_Phi_ref"].resize((hf["eul_Phi_ref"].shape[0]+1), axis=0)
			hf["eul_Phi_ref"][-1]=eul_Phi_ref
			hf["eul_Psi_ref"].resize((hf["eul_Psi_ref"].shape[0]+1), axis=0)
			hf["eul_Psi_ref"][-1]=eul_Psi_ref
			hf["Phi_T_ref"].resize((hf["Phi_T_ref"].shape[0]+1), axis=0)
			hf["Phi_T_ref"][-1]=Phi_T_ref

			# Append the surface response
			hf["xi_r_ref"].resize((hf["xi_r_ref"].shape[0]+1), axis=0)
			hf["xi_r_ref"][-1]=xi_r_ref
			hf["lag_L_ref"].resize((hf["lag_L_ref"].shape[0]+1), axis=0)
			hf["lag_L_ref"][-1]=lag_L_ref

	else:
		# terminate if tidal response history doesn't exist when it should
		assert pind==1, [print("Couldn't find tidal response history at pind={}.".format(pind)), sys.exit(1)]

		# creating datasets with maxshape=None allows them to be resized later
		with h5py.File("{}/tidal_response_history.h5".format(orbev_fodir), "a") as hf:
			# save the unique mode descriptors only once
			hf.create_dataset("l", data=lun)
			hf.create_dataset("m", data=mun)
			hf.create_dataset("k", data=kun)

			# save current orbital parameters
			hf.create_dataset("q", data=np.array([np.sum(q)])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			hf.create_dataset("R_a", data=np.array([np.sum(R_a)])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			hf.create_dataset("Omega_orb", data=np.array([np.sum(Omega_orb)])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			hf.create_dataset("e", data=np.array([np.sum(e)])[...,np.newaxis], chunks=True, maxshape=(None, 1))

			# save orbital evolution rates
			# change in argument of periastron
			hf.create_dataset("o_dot_freq", data=o_dots, chunks=True, maxshape=(None, *o_dots.shape[1:]))
			hf.create_dataset("o_dot", data=np.array([np.sum(o_dots)])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			# change in semi-major axis
			hf.create_dataset("a_dot_freq", data=a_dots, chunks=True, maxshape=(None, *a_dots.shape[1:]))
			hf.create_dataset("a_dot", data=np.array([np.sum(a_dots)])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			# change in eccentricity
			hf.create_dataset("e_dot_freq", data=e_dots, chunks=True, maxshape=(None, *e_dots.shape[1:]))
			hf.create_dataset("e_dot", data=np.array([np.sum(e_dots)])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			# change in rotational angular momentum
			hf.create_dataset("J_dot_freq", data=J_dots, chunks=True, maxshape=(None, *J_dots.shape[1:]))
			hf.create_dataset("J_dot", data=np.array([np.sum(J_dots)])[...,np.newaxis], chunks=True, maxshape=(None, 1))

			# save the fundamental parameters used to calculate the orbital evolution rates
			hf.create_dataset("Gbar_1", data=Gbar_1, chunks=True, maxshape=(None, *Gbar_1.shape[1:]))
			hf.create_dataset("Gbar_2", data=Gbar_2, chunks=True, maxshape=(None, *Gbar_2.shape[1:]))
			hf.create_dataset("Gbar_3", data=Gbar_3, chunks=True, maxshape=(None, *Gbar_3.shape[1:]))
			hf.create_dataset("Gbar_4", data=Gbar_4, chunks=True, maxshape=(None, *Gbar_4.shape[1:]))
			hf.create_dataset("cbar", data=cbar, chunks=True, maxshape=(None, *cbar.shape[1:]))
			hf.create_dataset("eul_Phi_ref", data=eul_Phi_ref, chunks=True, maxshape=(None, *eul_Phi_ref.shape[1:]))
			hf.create_dataset("eul_Psi_ref", data=eul_Psi_ref, chunks=True, maxshape=(None, *eul_Psi_ref.shape[1:]))
			hf.create_dataset("Phi_T_ref", data=Phi_T_ref, chunks=True, maxshape=(None, *Phi_T_ref.shape[1:]))

			# save the surface response (displacement and flux perturbation)
			hf.create_dataset("xi_r_ref", data=xi_r_ref, chunks=True, maxshape=(None, *xi_r_ref.shape[1:]))
			hf.create_dataset("lag_L_ref", data=lag_L_ref, chunks=True, maxshape=(None, *lag_L_ref.shape[1:]))

def get_kappa(m,k):
	if k==0:
		if m==0:
			kappa=0.5
		elif m>=1:
			kappa=1.0
		else:
			kappa=0.0
	elif k>0:
		kappa=1.0
	else:
		kappa=0.0

	return kappa

if __name__=="__main__":
	main()
