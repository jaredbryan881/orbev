import numpy as np
import pandas as pd
import h5py

import matplotlib.pyplot as plt
import matplotlib.cm as cm

def main():
	#base_diname="/home/jared/MIT/astero/gyre_HATP2/orbev/data/figure/live_planet"
	base_diname="/home/jaredbryan/MIT/astero/gyre_HATP2/orbev/data/long_time"
	base_star=["M1.36_Z0.02262", "M1.48_Z0.02262", "M1.24_Z0.02262", "M1.36_Z0.02996", "M1.36_Z0.01708"]

	suffix="loworot"

	fig,axs=plt.subplots(6,2,sharex=True,sharey=True,figsize=(5,10))

	#bins=np.linspace(2,9.5,40)
	bins=np.linspace(-14,-5,20)

	# load base trajectory
	try:
		print("{}/{}_orb0_{}_fwd/orbital_history.data".format(base_diname, base_star[0], suffix))
		base_h_fwd=pd.read_csv("{}/{}_orb0_{}_fwd/orbital_history.data".format(base_diname, base_star[0], suffix))
		base_h_rev=pd.read_csv("{}/{}_orb0_{}_rev/orbital_history.data".format(base_diname, base_star[0], suffix))

		base_h=pd.concat([base_h_rev[::-1], base_h_fwd])

		base_edot, base_adot = calculate_orbev(base_h)
		#base_Te,base_Ta=calculate_T_deviation(base_h)

		add_hist(axs[0], base_edot, base_adot, c='k', bins=bins)
		add_hist(axs[1], base_edot, base_adot, c='k', bins=bins)
		add_hist(axs[2], base_edot, base_adot, c='k', bins=bins)
		add_hist(axs[3], base_edot, base_adot, c='k', bins=bins)
		add_hist(axs[4], base_edot, base_adot, c='k', bins=bins)
		add_hist(axs[5], base_edot, base_adot, c='k', bins=bins)

	except Exception as ex:
		print(ex)
		print("Failed to load {}_orb{}".format(base_star, i))

	axs[0,0].set_yscale("log")

	# variation in eccentricity
	for i in range(1,4):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[0], i, suffix))
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[0], i, suffix))

			e_h=pd.concat([h_rev[::-1], h_fwd])

			#e_Te,a_Ta=calculate_T_deviation(e_h)
			#add_hist(axs[0], e_Te, a_Ta, c=cm.inferno((i+1)/5), bins=bins)

			e_edot, e_adot = calculate_orbev(e_h)
			add_hist(axs[0], e_edot, e_adot, c=cm.inferno((i+1)/5), bins=bins)

		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))
			continue

	# variation in semi-major axis
	for i in range(4,7):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[0], i, suffix))
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[0], i, suffix))

			a_h=pd.concat([h_rev[::-1], h_fwd])

			#a_Te,a_Ta=calculate_T_deviation(a_h)
			#add_hist(axs[1], e_Te, a_Ta, c=cm.inferno((i-2)/5))

			a_edot, a_adot = calculate_orbev(a_h)
			add_hist(axs[1], a_edot, a_adot, c=cm.inferno((i-2)/5), bins=bins)


		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))
			continue

	# variation in rotation rate
	for i in range(7,10):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[0], i, suffix))
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[0], i, suffix))

			Orot_h=pd.concat([h_rev[::-1], h_fwd])

			#Orot_Te,Orot_Ta=calculate_T_deviation(Orot_h)
			#add_hist(axs[2], Orot_Te, Orot_Ta, c=cm.inferno((i-5)/5))

			Orot_edot, Orot_adot = calculate_orbev(Orot_h)
			add_hist(axs[2], Orot_edot, Orot_adot, c=cm.inferno((i-5)/5), bins=bins)

		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))
			continue

	# variation in initial time
	for i in range(10,12):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[0], i, suffix))
		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))

		try:
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[0], i, suffix))
		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))

		if i==10:
			t_h=pd.concat([h_rev[::-1], h_fwd])
		elif i==11:
			t_h=h_rev[::-1]

		#t_Te,t_Ta=calculate_T_deviation(t_h)
		#add_hist(axs[3], t_Te, t_Ta, c=cm.inferno((i-8)/5))

		t_edot, t_adot = calculate_orbev(t_h)
		add_hist(axs[3], t_edot, t_adot, c=cm.inferno((i-8)/5), bins=bins)

	# variation in mass
	star_ind=1
	for i in range(12,14):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[star_ind], i, suffix))
			fwd_diff=(h_fwd.a.values[0]-base_h_fwd.a.values[0])
		except Exception as ex:
			print(ex)
			print("Failed to load {}_orb{}".format(base_star[star_ind], i))

		try:
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[star_ind], i, suffix))
			rev_diff=(h_rev.a.values[0]-base_h_rev.a.values[0])
		except Exception as ex:
			print(ex)
			print("Failed to load {}_orb{}".format(base_star[star_ind], i))
		star_ind+=1

		if i==12:
			M_h=pd.concat([h_rev[::-1], h_fwd])
		elif i==13:
			M_h=h_rev[::-1]

		#M_Te,M_Ta=calculate_T_deviation(M_h)
		#add_hist(axs[4], M_Te, M_Ta, c=cm.inferno((i-10)/5))

		M_edot, M_adot = calculate_orbev(M_h)
		add_hist(axs[4], M_edot, M_adot, c=cm.inferno((i-10)/5), bins=bins)

	# variation in metallicity
	for i in range(14,16):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[star_ind], i, suffix))
		except Exception as ex:
			print(ex)
			print("Failed to load {}_orb{}".format(base_star[star_ind], i))

		try:
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[star_ind], i, suffix))
		except Exception as ex:
			print(ex)
			print("Failed to load {}_orb{}".format(base_star[star_ind], i))
			continue
		star_ind+=1

		Z_h=pd.concat([h_rev[::-1], h_fwd])
		
		#Z_Te,Z_Ta=calculate_T_deviation(Z_h)
		#add_hist(axs[5], Z_Te, Z_Ta, c=cm.inferno((i-12)/5))

		Z_edot, Z_adot = calculate_orbev(Z_h)
		add_hist(axs[5], Z_edot, Z_adot, c=cm.inferno((i-12)/5), bins=bins)

	axs[0,0].set_ylabel("Density")
	axs[1,0].set_ylabel("Density")
	axs[2,0].set_ylabel("Density")
	axs[3,0].set_ylabel("Density")
	axs[4,0].set_ylabel("Density")
	axs[5,0].set_ylabel("Density")
	#axs[-1,0].set_xlabel(r"$T_e$")
	#axs[-1,1].set_xlabel(r"$T_a$")
	axs[-1,0].set_xlabel(r"$\dot{e}$")
	axs[-1,1].set_xlabel(r"$\dot{a}$")
	plt.tight_layout()
	#plt.savefig("Te_Ta_orbev.pdf")
	plt.savefig("edot_adot_orbev.pdf")
	plt.show()

def first_gt(arr, val):
	return np.argmax(arr>val)

def calculate_T_deviation(h, e_sig=0.00082, a_sig=0.00065):
	times=h.time.values
	e=h.e.values
	a=h.a.values

	Te=np.zeros(len(times))
	Ta=np.zeros(len(times))
	for (i_time, time) in enumerate(times):
		de=np.abs(e-e[i_time])[i_time:]
		Te_ind=first_gt(de, e_sig)
		Te[i_time]=times[i_time:][Te_ind]-time

		da=np.abs(a-a[i_time])[i_time:]
		Ta_ind=first_gt(da, a_sig)
		Ta[i_time]=times[i_time:][Ta_ind]-time

	return Te, Ta

def calculate_orbev(h):
	times=h.time.values
	e=h.e.values
	a=h.a.values

	de=e[1:]-e[:-1]
	da=a[1:]-a[:-1]
	dt=times[1:]-times[:-1]

	edot=de[dt!=0]/dt[dt!=0]
	adot=da[dt!=0]/dt[dt!=0]

	return edot, adot

def add_hist(axs, data1, data2, c, bins):
	axs[0].hist(np.log10(data1[np.nonzero(data1)]), facecolor=c, alpha=0.5, bins=bins, density=True)
	axs[1].hist(np.log10(data2[np.nonzero(data2)]), facecolor=c, alpha=0.5, bins=bins, density=True)

if __name__=="__main__":
	main()