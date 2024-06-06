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
	base_finame="."

	cur_star_save=""
	for fidir in os.listdir(base_finame)[9:]:
		print(fidir)
		cur_star_split=fidir.split("_")
		cur_star="{}_{}".format(cur_star_split[0], cur_star_split[1])

		if cur_star!=cur_star_save:
			print("Changing cur_star from {} to {}".format(cur_star_save, cur_star))
			cur_star_save="{}_{}".format(cur_star_split[0], cur_star_split[1])

			try:
				h=mr.MesaData("{}/history_{}.data".format(base_finame, cur_star))
			except:
				print("Failed to load {}".format(cur_star))
				continue

		try:
			calc_freq_scale(fidir, h)
		except Exception as ex:
			print(ex)
			print("Failed to calculate for {}".format(fidir))
			continue

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

def calc_freq_scale(fidir, h):
	ages=h.star_age
	R=h.radius*Rsun
	M_cur=float(fidir.split("_")[0][1:])
	M=M_cur*Msun
	Teff=10**(h.log_Teff)
	logg=h.log_g

	# calculate frequency scale at MESA timesteps
	fs=freq_scale(M,R)

	# load orbital history file
	h_orb=pd.read_csv("{}/orbital_history.data".format(fidir))
	time_orb=h_orb["time"]

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

	np.savetxt("{}/fs.txt".format(fidir), fs_orb)
	np.savetxt("{}/Ts.txt".format(fidir), T_orb)
	np.savetxt("{}/Rs.txt".format(fidir), R_orb)
	np.savetxt("{}/loggs.txt".format(fidir), logg_orb)

if __name__=="__main__":
	main()
