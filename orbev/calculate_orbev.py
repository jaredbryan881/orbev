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

	eul_Psi_ref_temp=np.zeros(len(eul_Psi_ref), dtype=np.complex128)
	for i in range(len(eul_Psi_ref)):
		eul_Psi_ref_temp[i] = eul_Psi_ref[i][0] + 1j*eul_Psi_ref[i][1]
	eul_phi=eul_Psi_ref_temp
	eul_phi.real-=Phi_T_ref

	n_l=len(np.unique(l))
	n_m=len(np.unique(m))
	n_k=len(np.unique(k))

	o_dots=np.zeros((n_l,n_m,n_k))
	a_dots=np.zeros((n_l,n_m,n_k))
	e_dots=np.zeros((n_l,n_m,n_k))
	J_dots=np.zeros((n_l,n_m,n_k))

	for (i1,l_cur) in enumerate(np.unique(l)):
		for (i2,m_cur) in enumerate(np.unique(m)):
			for (i3,k_cur) in enumerate(np.unique(k)):
				ind = i2*n_k + i3

				c=cbar[ind]
				kappa=get_kappa(m_cur,k_cur)
				if (c==0) or (kappa==0):
					continue

				F = -(1/2)*(np.sqrt(4*np.pi)*eul_phi[ind]/(eps_tide*c) + 1)
				gamma = np.arctan2(np.imag(F), np.real(F))

				o_dots[i1,i2,i3] = 4*Omega_orb*q*(R_a)**(l_cur+3)*kappa*np.abs(F)*np.cos(gamma)*Gbar_1[ind]
				a_dots[i1,i2,i3] = 4*Omega_orb*(q/R_a)*(R_a)**(l_cur+3)*kappa*np.abs(F)*np.sin(gamma)*Gbar_2[ind]
				e_dots[i1,i2,i3] = 4*Omega_orb*q*(R_a)**(l_cur+3)*kappa*np.abs(F)*np.sin(gamma)*Gbar_3[ind]
				J_dots[i1,i2,i3] = 4*Omega_orb*q**2/np.sqrt(R_a*(1+q))*(R_a)**(l_cur+3)*kappa*np.abs(F)*np.sin(gamma)*Gbar_4[ind]

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