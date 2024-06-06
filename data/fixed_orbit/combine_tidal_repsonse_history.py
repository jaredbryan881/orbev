import numpy as np
import h5py
import matplotlib.pyplot as plt 

# constants
G = 6.67430e-11  # [m3 kg-1 s-2]
Rsun = 695700000 # [m]
Msun = 1.989e30  # [kg]
au = 1.496e+11   # [m]

def main():
	# stellar parameters
	Ms=[1.36, 1.36136, 1.3634, 1.3668, 1.3702, 1.3736, 1.394, 1.428, 1.36, 1.36, 1.36, 1.36, 1.36, 1.36, 1.36]
	Zs=[0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02264262, 0.02267655, 0.0227331, 0.02278965, 0.0228462, 0.0231855, 0.023751]

	# group stellar parameters over orbital configurations
	n_times_max=7096
	for i_star in range(15):
		cur_star="M{}_Z{}".format(Ms[i_star], Zs[i_star])
		# load tidal response rates
		finame="tidal_response_history_{}.h5".format(i_star)
		with h5py.File(finame, "r") as hf:
			e_dot_freq=hf["e_dot_freq"][:]
			a_dot_freq=hf["a_dot_freq"][:]
			e=hf["e"][1]
			Omega_orb=hf["Omega_orb"][1]
		a = OmegaOrb_to_a(Omega_orb, Ms[i_star])

		n_times=e_dot_freq.shape[0]

		# load dimensionalization parameters
		fs=np.loadtxt("fs_{}.txt".format(cur_star))
		Rs=np.loadtxt("Rs_{}.txt".format(cur_star))

		# initialize save arrays
		if i_star==0:
			e_dot_freq_allstar=np.zeros((15, n_times_max, 2, 51))
			a_dot_freq_allstar=np.zeros((15, n_times_max, 2, 51))
			e_allstar=np.zeros(15)
			a_allstar=np.zeros(15)
			Omega_orb_allstar=np.zeros(15)

		e_dot_freq_allstar[i_star, :n_times, :, :]=(e_dot_freq.T*fs[:n_times]).T
		a_dot_freq_allstar[i_star, :n_times, :, :]=(a_dot_freq.T*fs[:n_times]*Rs[:n_times]/au).T
		e_allstar[i_star]=e
		a_allstar[i_star]=a
		Omega_orb_allstar[i_star]=Omega_orb

		e_dot_freq_allstar[i_star, n_times:, :, :]=np.nan
		a_dot_freq_allstar[i_star, n_times:, :, :]=np.nan

	with h5py.File("tidal_response_history_allstar.h5", "w") as hf:
		hf.create_dataset("e_dot_freq", data=e_dot_freq_allstar)
		hf.create_dataset("a_dot_freq", data=a_dot_freq_allstar)
		hf.create_dataset("e", data=e_allstar)
		hf.create_dataset("a", data=a_allstar)
		hf.create_dataset("Omega_orb", data=Omega_orb_allstar)

def OmegaOrb_to_a(OmegaOrb, M):
	"""Convert orbital frequency to semi-major axis

	Arguments
	---------
	:param OmegaOrb: float
		Orbital frequency [cyc/day]

	Returns
	-------
	:return a: float
		Semi-major axis [au]
	"""
	T=1/OmegaOrb # [day]
	T*=86400 # [s]
	a=((G*M*Msun*T**2)/(4*np.pi**2))**(1/3) # [m]
	a/=au # [au]

	return a

if __name__=="__main__":
	main()