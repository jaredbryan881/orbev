import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

def main():
	# load command line arguments
	pind=int(sys.argv[1])
	cur_dir=sys.argv[2]

	with h5py.File("./profile{}/tide_orbit.h5".format(pind-1), "r") as hf:
		l=hf["l"][:]
		m=hf["m"][:]
		k=hf["k"][:]

		q=hf["q"][:]
		R_a=hf["R_a"][:]
		Omega_orb=hf["Omega_orb"][:]
		e=hf["e"][:]

		# secular evolution coefficients (G factors)
		Gbar_1=hf["Gbar_1"][:]
		Gbar_2=hf["Gbar_2"][:]
		Gbar_3=hf["Gbar_3"][:]
		Gbar_4=hf["Gbar_4"][:]

		cbar=hf["cbar"][:]

		eul_Psi_ref=hf["eul_Psi_ref"][:]
		Phi_T_ref=hf["Phi_T_ref"][:]

	Omega_orb=Omega_orb[0]
	R_a=R_a[0]
	q=q[0]
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

	o_dots=np.zeros((n_l,n_m,n_k))
	a_dots=np.zeros((n_l,n_m,n_k))
	e_dots=np.zeros((n_l,n_m,n_k))
	J_dots=np.zeros((n_l,n_m,n_k))

	print(len(l))
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

	with h5py.File("./profile{}/tide_orbit.h5".format(pind-1), "a") as hf:
		hf.create_dataset("o_dot_freq", data=o_dots)
		hf.create_dataset("o_dot", data=[np.sum(o_dots)])

		hf.create_dataset("a_dot_freq", data=a_dots)
		hf.create_dataset("a_dot", data=[np.sum(a_dots)])

		hf.create_dataset("e_dot_freq", data=e_dots)
		hf.create_dataset("e_dot", data=[np.sum(e_dots)])

		hf.create_dataset("J_dot_freq", data=J_dots)
		hf.create_dataset("J_dot", data=[np.sum(J_dots)])

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