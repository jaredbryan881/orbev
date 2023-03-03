import numpy as np
import h5py
import sys
import os

import pygyre as pg

def main():
	# load command line arguments
	orbev_fodir=sys.argv[1]
	pind=int(sys.argv[2])

	s=pg.read_output("{}/tide_orbit.h5".format(orbev_fodir))

	# Extract the first set of responses
	sg = s.group_by('id').groups[0]

	Omega_orb = sg['Omega_orb']  
	R_a = sg['R_a']
	q = sg['q']
	e = sg['e']

	eps_T = (R_a**3)*q

	l = sg['l']
	m = sg['m']
	k = sg['k']

	cbar = sg['cbar']

	Fbar = -(1/2)*sg['eul_Phi_ref']/(cbar*eps_T)

	x = sg['x_ref']

	Gbar_1 = sg['Gbar_1']
	Gbar_2 = sg['Gbar_2']
	Gbar_3 = sg['Gbar_3']
	Gbar_4 = sg['Gbar_4']

	# radial displacement perturbation at reference location
	xi_r_ref=sg["xi_r_ref"]
	xi_r_ref=xi_r_ref.real + 1j*xi_r_ref.imag
	# Lagrangian radiative luminosity perturbation at reference location
	lag_L_ref=sg["lag_L_ref"]
	lag_L_ref=lag_L_ref.real + 1j*lag_L_ref.imag

	kap = np.empty(len(l))
	for i in range(len(kap)):
		kap[i]=get_kappa(m[i], k[i])

	# Argument of periastron (units of radians per dynamical timescale)
	o_dots = 4*Omega_orb*q * (R_a)**(l+3) * (x)**(l+1) * kap * Fbar.real * Gbar_1
	# Semi-major axis (units of R per dynamical timescale)
	a_dots = 4*Omega_orb*(q/R_a) * (R_a)**(l+3) * (x)**(l+1) * kap * Fbar.imag * Gbar_2
	# Eccentricity (units of per dynamical timescale)
	e_dots = 4*Omega_orb*q * (R_a)**(l+3) * (x)**(l+1) * kap * Fbar.imag * Gbar_3
	# Angular momentum (units of GM^2/R)
	J_dots = 4*(q**2)*R_a * (R_a)**(l+3) * (x)**(l+1) * kap * Fbar.imag * Gbar_4

	# reshape arrays to keep only m=-2,0,2
	# odots
	o_dots=np.reshape(o_dots, (1,5,51))
	o_dots=o_dots[:,[0,2,4],:]
	# adots
	a_dots=np.reshape(a_dots, (1,5,51))
	a_dots=a_dots[:,[0,2,4],:]
	# edots
	e_dots=np.reshape(e_dots, (1,5,51))
	e_dots=e_dots[:,[0,2,4],:]
	# Jdots
	J_dots=np.reshape(J_dots, (1,5,51))
	J_dots=J_dots[:,[0,2,4],:]

	# 
	cbar=np.reshape(cbar, (1,5,51))
	cbar=cbar[:,[0,2,4],:]
	#
	Fbar=np.reshape(Fbar, (1,5,51))
	Fbar=Fbar[:,[0,2,4],:]
	#
	Gbar_1=np.reshape(Gbar_1, (1,5,51))
	Gbar_1=Gbar_1[:,[0,2,4],:]
	Gbar_2=np.reshape(Gbar_2, (1,5,51))
	Gbar_2=Gbar_2[:,[0,2,4],:]
	Gbar_3=np.reshape(Gbar_3, (1,5,51))
	Gbar_3=Gbar_3[:,[0,2,4],:]
	Gbar_4=np.reshape(Gbar_4, (1,5,51))
	Gbar_4=Gbar_4[:,[0,2,4],:]

	# radial displacement
	xi_r_ref=np.reshape(xi_r_ref, (1,5,51))
	xi_r_ref=xi_r_ref[:,[0,2,4],:]
	# lagrangian radiative luminosity perturbation
	lag_L_ref=np.reshape(lag_L_ref, (1,5,51))
	lag_L_ref=lag_L_ref[:,[0,2,4],:]

	# Write data
	if os.path.exists("{}/tidal_response_history.h5".format(orbev_fodir)):
		# Resize datasets and append current values
		with h5py.File("{}/tidal_response_history.h5".format(orbev_fodir), "a") as hf:
			# Append current orbital parameters
			# ratio of secondary mass to primary mass
			hf["q"].resize((hf["q"].shape[0]+1), axis=0)
			hf["q"][-1]=np.array([q[0]])
			# ratio of primary radius to semi-major axis
			hf["R_a"].resize((hf["R_a"].shape[0]+1), axis=0)
			hf["R_a"][-1]=np.array([R_a[0]])
			# orbital frequency
			hf["Omega_orb"].resize((hf["Omega_orb"].shape[0]+1), axis=0)
			hf["Omega_orb"][-1]=np.array([Omega_orb[0]])
			# eccentricity
			hf["e"].resize((hf["e"].shape[0]+1), axis=0)
			hf["e"][-1]=np.array([e[0]])

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
			hf["Fbar"].resize((hf["Fbar"].shape[0]+1), axis=0)
			hf["Fbar"][-1]=Fbar

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
			hf.create_dataset("l", data=np.unique(l))
			hf.create_dataset("m", data=np.unique(m))
			hf.create_dataset("k", data=np.unique(k))

			# save current orbital parameters
			hf.create_dataset("q", data=np.array([np.sum(q[0])])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			hf.create_dataset("R_a", data=np.array([np.sum(R_a[0])])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			hf.create_dataset("Omega_orb", data=np.array([np.sum(Omega_orb[0])])[...,np.newaxis], chunks=True, maxshape=(None, 1))
			hf.create_dataset("e", data=np.array([np.sum(e[0])])[...,np.newaxis], chunks=True, maxshape=(None, 1))

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
			hf.create_dataset("Fbar", data=Fbar, chunks=True, maxshape=(None, *Fbar.shape[1:]))

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
