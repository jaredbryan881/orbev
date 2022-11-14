import numpy as np

# constants
G = 6.67430e-11  # [m3 kg-1 s-2]
Rsun = 695700000 # [m]
Msun = 1.989e30  # [kg]
au = 1.496e+11   # [m]

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

def a_to_OmegaOrb(a, M):
	"""Convert semi-major axis to orbital frequency

	Arguments
	---------
	:param a: float
		Semi-major axis [au]

	Returns
	-------
	:param OmegaOrb: float
		Orbital frequency [cyc/day]
	"""
	a*=au # [m]
	T=2*np.pi*np.sqrt(a**3/(G*M*Msun)) # [s]
	T/=86400 # [day]
	OmegaOrb=1/T # [cyc/day]

	return OmegaOrb
