import numpy as np

# Define constants
Rsun=695700 # km
Msun=1.9891e30 # kg
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
	R=1.7039200 # Rsun
	R_sig=[0.0528432, 0.0804382] # Rsun

	# Stellar mass: TICv8
	M=1.29 # Msun
	M_sig=[0.158558, 0.224551] # Msun

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
	Mp_sig=[0.20, 0.19] # Mjup

	Orot = Orot_from_vR(v, R) # cyc/day
	Oorb = 1/Porb # cyc/day

	#########################
	# Covary all parameters #
	#########################
	# get rotational frequency
	vs = v_sig*np.random.randn(N) + v
	Rs = np.mean(R_sig)*np.random.randn(N) + R
	Orots = Orot_from_vR(vs, Rs) # [cyc/day]
	Orot_sig = np.std(Orots)

	# get orbital frequency
	Porbs = np.mean(Porb_sig)*np.random.randn(N) + Porb
	Oorbs = 1/Porbs
	Oorb_sig=np.std(Oorbs)

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
	q=np.mean(qs)
	q_sig=np.std(qs)

	# get Zs
	MHs = MH_sig*np.random.randn(N) + MH
	Zs = MH2Z(MHs)
	Z = np.mean(Zs)

	#np.savetxt("params.txt", [es, Oorbs, Orots, ages, qs, Ms, Zs])
	print("e={:3f},     min={:3f}, max={:3f}".format(np.mean(es), np.min(es), np.max(es)))
	print("Oorbs={:3f}, min={:3f}, max={:3f}".format(np.mean(Oorbs), np.min(Oorbs), np.max(Oorbs)))
	print("Orots={:3f}, min={:3f}, max={:3f}".format(np.mean(Orots), np.min(Orots), np.max(Orots)))
	print("ages={:3f},  min={:3f}, max={:3f}".format(np.mean(ages), np.min(ages), np.max(ages)))
	print("qs={:3f},    min={:3f}, max={:3f}".format(np.mean(qs), np.min(qs), np.max(qs)))
	print("Ms={:3f},    min={:3f}, max={:3f}".format(np.mean(Ms)/Msun, np.min(Ms)/Msun, np.max(Ms)/Msun))
	print("Zs={:3f},    min={:3f}, max={:3f}".format(np.mean(Zs), np.min(Zs), np.max(Zs)))

	np.savetxt("params_low_orot_ensemble.txt", [es, Oorbs, Orots, ages, Ms, Zs])

	Orots=Oorbs*psOrot(es)
	np.savetxt("params_ps_orot_ensemble.txt", [es, Oorbs, Orots, ages, Ms, Zs])

	##########################
	# Vary Single parameters #
	##########################
	es_vary    = [e+e_sig[1], e-2*e_sig[0], e+3*e_sig[1]] # 3
	Oorbs_vary = [Oorb+Oorb_sig, Oorb-2*Oorb_sig, Oorb+3*Oorb_sig] # 3
	Orots_vary = [Orot+Orot_sig, Orot-2*Orot_sig, Orot+3*Orot_sig] # 3
	ages_vary  = [age+age_sig, age-age_sig] # 2
	Ms_vary    = [M+M_sig[1], M-M_sig[0]] # 2
	Zs_vary    = [0.02*(10**(MH+MH_sig)), 0.02*(10**(MH-MH_sig))] # 3

	es_save    = [e]        + [e for i in range(0)]         + es_vary    + [e for i in range(12)]
	Oorbs_save = [Oorb]     + [Oorb for i in range(3)]      + Oorbs_vary + [Oorb for i in range(9)]
	Orots_save = [Orot]     + [Orot for i in range(6)]      + Orots_vary + [Orot for i in range(6)]
	ages_save  = [age]      + [age for i in range(9)]       + ages_vary  + [age for i in range(4)]
	Ms_save    = [M]        + [M for i in range(11)]        + Ms_vary    + [M for i in range(2)]
	Zs_save    = [MH2Z(MH)] + [MH2Z(MH) for i in range(13)] + Zs_vary    + [MH2Z(MH) for i in range(0)]

	# HAT-P-2 Orot
	Orots_vary = [Orot+Orot_sig, Orot-2*Orot_sig, Orot+3*Orot_sig]
	Orots_save = [Orot]     + [Orot for i in range(6)]      + Orots_vary + [Orot for i in range(6)]
	np.savetxt("params_low_orot.txt", [es_save, Oorbs_save, Orots_save, ages_save, Ms_save, Zs_save])

	# Pseudosynchronous Orot
	Orot = Oorb*psOrot(e)
	Orots_vary = [Orot+Orot_sig, Orot-2*Orot_sig, Orot+3*Orot_sig]
	Orots_save = [Orot]     + [Orot for i in range(6)]      + Orots_vary + [Orot for i in range(6)]
	np.savetxt("params_ps_orot.txt", [es_save, Oorbs_save, Orots_save, ages_save, Ms_save, Zs_save])

def Orot_from_vR(v, R):
	Prot = (2*np.pi*R*Rsun)/v # s
	Prot/=86400 # day
	Orot = 1/Prot # cyc/day
	return Orot

def psOrot(e):
	num = (1 + (15/2)*(e**2) + (45/8)*(e**4) + (5/16)*(e**6))
	denom = (1 + 3*(e**2) + (3/8)*(e**4))*((1-(e**2))**(3/2))
	return num/denom

def MH2Z(MH):
	return 0.02*(10**MH)

if __name__=="__main__":
	main()