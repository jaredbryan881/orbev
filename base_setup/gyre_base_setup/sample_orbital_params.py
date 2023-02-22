import numpy as np

# Define constants
Rsun=695700 # km
Msun=1.9891e30  # kg
Mjup=1.899e27 # kg

def main():
	N=1000

	######################
	# Stellar parameters #
	######################
	# Surface velocity: TICv8
	# assume vsini=v
	v=20.8 # km/s
	v_sig=0.3 # km/s

	# Stellar radius: TICv8
	R=1.7039200
	R_sig=[0.0528432, 0.0804382]

	# Stellar mass: TICv8
	M=1.29
	M_sig=[0.158558, 0.224551]

	# Stellar metallicity: TICv8
	MH = 0.0535135 # dex
	MH_sig = 0.0406867 # dex

	# Stellar age: Bonomo et al. (2017)
	age=2.60 # Gyr
	age_sig=0.5 # Gyr

	######################
	# Orbital parameters #
	######################
	# Orbital period: Bonomo et al. (2017)
	Porb=5.6334754 # day
	Porb_sig=0.0000026 # day

	# Orbital eccentricity: Bonomo et al. (2017)
	e=0.50833
	e_sig=[0.00075, 0.00082]

	# Stellar mass: Bonomo et al. (2017)
	Mp=8.70 # Mjup
	Mp_sig=[0.20, 0.19]

	Orot = Orot_from_vR(v, R)
	Oorb = 1/Porb

	#####################
	# Sample parameters #
	#####################
	# get rotational frequency
	vs = v_sig*np.random.randn(N) + v
	Rs = np.mean(R_sig)*np.random.randn(N) + R
	Orots = Orot_from_vR(vs, Rs) # [cyc/day]

	# get orbital frequency
	Porbs = np.mean(Porb_sig)*np.random.randn(N) + Porb
	Oorbs = 1/Porbs

	# get orbital eccentricity
	es = np.mean(e_sig)*np.random.randn(N) + e

	# get stellar age
	ages = age_sig*np.random.randn(N) + age

	# get q: Mp/Ms
	Ms = np.mean(M_sig)*np.random.randn(N) + M
	Ms*=Msun
	Mps = np.mean(Mp_sig)*np.random.randn(N) + Mp
	Mps*=Mjup
	qs = Mps/Ms

	# get Zs
	MHs = MH_sig*np.random.randn(N) + MH
	Zs = 0.02*10**MHs

	np.savetxt("params.txt", [es, Oorbs, Orots, ages, qs, Ms, Zs])
	print("e={:3f},     min={:3f}, max={:3f}".format(np.mean(es), np.min(es), np.max(es)))
	print("Oorbs={:3f}, min={:3f}, max={:3f}".format(np.mean(Oorbs), np.min(Oorbs), np.max(Oorbs)))
	print("Orots={:3f}, min={:3f}, max={:3f}".format(np.mean(Orots), np.min(Orots), np.max(Orots)))
	print("ages={:3f},  min={:3f}, max={:3f}".format(np.mean(ages), np.min(ages), np.max(ages)))
	print("qs={:3f},    min={:3f}, max={:3f}".format(np.mean(qs), np.min(qs), np.max(qs)))
	print("Zs={:3f},    min={:3f}, max={:3f}".format(np.mean(Zs), np.min(Zs), np.max(Zs)))

def Orot_from_vR(v, R):
	Prot = (2*np.pi*R*Rsun)/v # s
	Prot/=86400 # day
	Orot = 1/Prot # cyc/day
	return Orot

if __name__=="__main__":
	main()