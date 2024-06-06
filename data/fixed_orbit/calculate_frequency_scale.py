import os
import numpy as np
import mesa_reader as mr
import pandas as pd
from scipy.interpolate import interp1d

# constants
G = 6.67430e-11  # [m3 kg-1 s-2]
Rsun = 695700000 # [m]
Msun = 1.989e30  # [kg]
au = 1.496e+11   # [m]

def main():
	Ms=[1.36, 1.36136, 1.3634, 1.3668, 1.3702, 1.3736, 1.394, 1.428, 1.36, 1.36, 1.36, 1.36, 1.36, 1.36, 1.36]
	Zs=[0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02264262, 0.02267655, 0.0227331, 0.02278965, 0.0228462, 0.0231855, 0.023751]
	for i, (M,Z) in enumerate(zip(Ms, Zs)):
		cur_star="M{}_Z{}".format(M,Z)
		print(cur_star)
		h_star=mr.MesaData("history_{}.data".format(cur_star))
		fs_orb, T_orb, R_orb, logg_orb = calc_freq_scale(h_star, M)

		np.savetxt("fs_{}.txt".format(cur_star), fs_orb)
		np.savetxt("Ts_{}.txt".format(cur_star), T_orb)
		np.savetxt("Rs_{}.txt".format(cur_star), R_orb)
		np.savetxt("loggs_{}.txt".format(cur_star), logg_orb)

def freq_scale(M, R):
	"""Calculate the freqency scale used for (non)dimensionalization.

	Arguments
	--------- 
	:param M: float
		Stellar mass [kg]
	:param R: float
		Stellar radius [m]
	"""
	# dimensionalizing constant
	freq_scale = np.sqrt(G*M/R**3) # rad/sec
	freq_scale*=(86400/(2*np.pi)) # cyc/day

	return freq_scale

def calc_freq_scale(h_star, M):
	ages=h_star.star_age
	R=h_star.radius*Rsun
	M*=Msun
	Teff=10**(h_star.log_Teff)
	logg=h_star.log_g

	# calculate frequency scale at MESA timesteps
	fs=freq_scale(M,R)

	# load orbital history file
	time_orb=np.linspace(2.6e9, 2.6e9 + 1e5*7609, 7609)

	# guard against overflow at the ends due to PMS and RGB
	time_orb[time_orb>ages[-1]]=ages[-1]
	time_orb[time_orb<ages[0]]=ages[0]

	# interpolate frequency scale to orbital history times
	f=interp1d(ages, fs)
	fs_orb=f(time_orb)

	# interpolate Teff to orbital history times
	f=interp1d(ages, Teff)
	T_orb=f(time_orb)

	# interpolate R to orbital history times
	f=interp1d(ages, R)
	R_orb=f(time_orb)

	# interpolate R to orbital history times
	f=interp1d(ages, logg)
	logg_orb=f(time_orb)

	return fs_orb, T_orb, R_orb, logg_orb

if __name__=="__main__":
	main()
