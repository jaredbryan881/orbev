import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

G = 6.67430e-11  # [m3 kg-1 s-2]
Rsun = 695700000 # [m]
Msun = 1.989e30  # [kg]
Rjup = 6.9911e7  # [m]
Mjup = 1.899e27  # [kg]
au = 1.496e+11   # [m]

def main():
	###
	### stellar parameters
	###
	M=1.29 # Msun
	M_sig=[0.158558, 0.224551] # Msun

	R=1.7039200 # Rsun
	R_sig=[0.0528432, 0.0804382] # Rsun

	###
	### planetary parameters
	###
	Mp=8.70 # Mjup
	Mp_sig=[0.20, 0.19] # Mjup

	Rp = 1.157 # Rjup
	Rp_sig = [0.063, 0.073]

	###
	### orbital parameters
	###
	Porb=5.6334754 # day
	Porb_sig=0.0000026 # day
	Omega_orb=1/Porb
	a0=OmegaOrb_to_a(Omega_orb, M)

	e0=0.50833
	e0_sig=[0.00075, 0.00082]

	###
	### simulation parameters
	###
	dt = 1e5 # yr
	
	t0=2.60 # Gyr
	t0_sig=0 #0.5 # Gyr

	time_reversed=False
	store_freq = int(1e2)
	n_steps = int(1e4)

	#########################
	# Covary all parameters #
	#########################
	N=10
	# get rotational frequency
	Rs = np.mean(R_sig)*np.random.randn(N) + R

	Rp = 1.157 # Rjup
	Rps = np.mean(Rp_sig)*np.random.randn(N) + Rp

	# get orbital frequency
	Porbs = np.mean(Porb_sig)*np.random.randn(N) + Porb
	Oorbs = 1/Porbs

	# get orbital eccentricity
	e0s = np.mean(e0_sig)*np.random.randn(N) + e0

	# get stellar age
	t0s = t0_sig*np.random.randn(N) + t0

	Ms = np.mean(M_sig)*np.random.randn(N) + M
	Mps = np.mean(Mp_sig)*np.random.randn(N) + Mp

	a0s = OmegaOrb_to_a(Oorbs, M)

	npts=int(n_steps/store_freq)

	tides='both'

	Qs = np.linspace(5,7,N)
	Qps = np.linspace(5,7,N)

	# vary each quantity individually
	quant_names=["Qp", "Q", "e0", "a0", "M", "R", "Mp", "Rp", "all"]
	for quant_name in quant_names:
		print(quant_name)
		t_arr=np.empty((N,2*npts))
		a_arr=np.empty((N,2*npts))
		e_arr=np.empty((N,2*npts))
		adot_arr=np.empty((N,2*npts))
		edot_arr=np.empty((N,2*npts))

		fig1,axs1=plt.subplots(2,1,sharex=True)
		for i in range(N):
			print(i)

			t0_samp = np.mean(t0s)
			e0_samp = np.mean(e0s)
			a0_samp = np.mean(a0s)
			M_samp  = np.mean(Ms)
			R_samp  = np.mean(Rs)
			Mp_samp = np.mean(Mps)
			Rp_samp = np.mean(Rps)
			Qp_samp = 10**6.5
			Q_samp  = 10**5.5
			if quant_name=="e0":
				e0_samp=e0s[i]
			elif quant_name=="a0":
				a0_samp=a0s[i]
			elif quant_name=="M":
				M_samp=Ms[i]
			elif quant_name=="R":
				R_samp=Rs[i]
			elif quant_name=="Mp":
				Mp_samp=Mps[i]
			elif quant_name=="Rp":
				Rp_samp=Rps[i]
			elif quant_name=="Qp":
				Qp_samp=10**Qps[i]
			elif quant_name=="Q":
				Q_samp=10**Qs[i]
			elif quant_name=="all":
				t0_samp = t0s[i]
				e0_samp = e0s[i]
				a0_samp = a0s[i]
				M_samp = Ms[i]
				R_samp = Rs[i]
				Mp_samp = Mps[i]
				Rp_samp = Rps[i]
				Qp_samp = 10**Qps[i]
				Q_samp = 10**Qs[i]

			color='grey'

			t_arr[i,:npts], a_arr[i,:npts], e_arr[i,:npts], adot_arr[i,:npts], edot_arr[i,:npts] = evolve_orbit(t0_samp, e0_samp, a0_samp, dt, M_samp, R_samp, Q_samp, Mp_samp, Rp_samp, Qp_samp, time_reversed=True, tides=tides, n_steps=n_steps)
			t_arr[i,npts:], a_arr[i,npts:], e_arr[i,npts:], adot_arr[i,npts:], edot_arr[i,npts:] = evolve_orbit(t0_samp, e0_samp, a0_samp, dt, M_samp, R_samp, Q_samp, Mp_samp, Rp_samp, Qp_samp, tides=tides, n_steps=n_steps)

			axs1[0].scatter(t0_samp/1e9, e0_samp, c='k')
			axs1[0].plot(t_arr[i]/1e9, e_arr[i], c=color, alpha=0.2)

			axs1[1].scatter(t0_samp/1e9, a0_samp, c='k')
			axs1[1].plot(t_arr[i]/1e9, a_arr[i], c=color, alpha=0.2)

		axs1[0].plot(np.nanmedian(t_arr, axis=0)/1e9, np.nanmedian(e_arr, axis=0), c='k')
		axs1[1].plot(np.nanmedian(t_arr, axis=0)/1e9, np.nanmedian(a_arr, axis=0), c='k')
		
		axs1[1].set_xlabel("Time [Gyr]", fontsize=14)
		axs1[0].set_ylabel("e", fontsize=14)
		axs1[1].set_ylabel("a [au]", fontsize=14)
		axs1[0].set_xlim(np.min(t_arr)/1e9, np.max(t_arr)/1e9)
		axs1[0].set_ylim(0, 0.75)
		axs1[1].set_ylim(0, 0.1)
		axs1[0].set_title(quant_name, fontsize=14)
		plt.tight_layout()
		#plt.savefig("figs/{}_tides_{}.pdf".format(tides, quant_name))
		#plt.savefig("figs/{}_tides_{}.png".format(tides, quant_name))
		#plt.close()
		plt.show()

def evolve_orbit(t0, e0, a0, dt, M, R, Q, Mp, Rp, Qp, n_steps=int(1e5), time_reversed=False, store_freq=int(1e2), tides='both'):
	a = []
	e = []
	t = []
	edots = []
	adots = []

	t_cur = t0
	e_cur = e0
	a_cur = a0
	for t_index in range(n_steps):
		# calculate change in eccentricity
		dedt = edot(a_cur*au, e_cur, Mp*Mjup, Rp*Rjup, Qp, M*Msun, R*Rsun, Q, tides=tides) # [1/s]
		# dimensionalize the value
		dedt = np.array(dedt)
		dedt = dedt*((86400*365)/(2*np.pi)) # [1/yr]

		# calculate change in semimajor axis
		dadt = adot(a_cur*au, e_cur, Mp*Mjup, Rp*Rjup, Qp, M*Msun, R*Rsun, Q, tides=tides) # [m/s]
		# dimensionalize the value
		dadt = np.array(dadt)
		dadt = dadt*((86400*365)/(2*np.pi)) # [m/yr]
		dadt = dadt/au # [au/yr]

		# update the values
		if time_reversed:
			e_cur -= dt*dedt.sum()
			a_cur -= dt*dadt.sum()
			t_cur -= dt
		else:
			e_cur += dt*dedt.sum()
			a_cur += dt*dadt.sum()
			t_cur += dt

		if t_index%store_freq == 0:
			# store the values
			t.append(t_cur)
			a.append(a_cur)
			e.append(e_cur)

			# store the rates
			edots.append(dedt.sum())
			adots.append(dadt.sum())

	if time_reversed:
		return t[::-1], a[::-1], e[::-1], adots[::-1], edots[::-1]
	else:
		return t, a, e, adots, edots

def edot(a, e, Mp, Rp, Qp, M, R, Q, tides='both'):
	# planetary tides
	p_tides = (63/4) * ((G*M**3)**(1/2)) * ((Rp**5)/(Qp*Mp))
	p_tides*=(-e * (a**(-13/2)))

	# stellar tides
	s_tides = (171/16) * ((G/M)**(1/2)) * ((Mp*R**5)/(Q))
	s_tides*=(-e * (a**(-13/2)))

	if tides=='planet':
		return p_tides, 0
	elif tides=='star':
		return 0, s_tides
	else:
		return p_tides, s_tides

def adot(a, e, Mp, Rp, Qp, M, R, Q, tides='both'):
	# planetary tides
	p_tides = (63/2) * ((G*M**3)**(1/2)) * ((Rp**5)/(Qp*Mp)) * (e**2)
	p_tides*=(-(a**(-11/2)))

	# stellar tides
	s_tides = (9/2) * ((G/M)**(1/2)) * ((Mp*R**5)/(Q))
	s_tides*=(-(a**(-11/2)))

	if tides=='planet':
		return p_tides, 0
	elif tides=='star':
		return 0, s_tides
	else:
		return p_tides, s_tides

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