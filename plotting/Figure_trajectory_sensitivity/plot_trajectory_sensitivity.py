import numpy as np
import pandas as pd
import h5py

import matplotlib.pyplot as plt

def main():
	#base_diname="/home/jared/MIT/astero/gyre_HATP2/orbev/data/figure/live_planet"
	base_diname="/home/jaredbryan/MIT/astero/gyre_HATP2/orbev/data/long_time"
	base_star=["M1.36_Z0.02262", "M1.48_Z0.02262", "M1.24_Z0.02262", "M1.36_Z0.02996", "M1.36_Z0.01708"]

	colors=['k'] + 3*['green'] + 3*['steelblue'] + 3*['crimson'] + 2*['purple'] + 2*['chocolate'] + 2*["slategray"]
	linestyles=['solid'] + 3*['solid', 'dashed', (0, (1,0.5))] + 3*['solid', "dashed"]
	linewidths=16*[2]

	suffix="highorot" # or loworot

	fig,axs=plt.subplots(6,2,sharex=True,sharey='col',figsize=(10,10))
	axs[-1,0].set_xlabel("Time [Myr]", fontsize=14)
	axs[-1,1].set_xlabel("Time [Myr]", fontsize=14)

	#axs[0,0].set_ylim(0.3, 0.7)
	#axs[0,1].set_ylim(0.05, 0.095)

	# load base trajectory
	try:
		print("{}/{}_orb0_{}_fwd/orbital_history.data".format(base_diname, base_star[0], suffix))
		base_h_fwd=pd.read_csv("{}/{}_orb0_{}_fwd/orbital_history.data".format(base_diname, base_star[0], suffix))
		base_h_rev=pd.read_csv("{}/{}_orb0_{}_rev/orbital_history.data".format(base_diname, base_star[0], suffix))
	except Exception as ex:
		print(ex)
		print("Failed to load {}_orb{}".format(base_star[i], i))

	# variation in eccentricity
	for i in range(1,4):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[0], i, suffix))
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[0], i, suffix))
		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))
			continue
		add_trajectory(axs[0], h_fwd, c=colors[i], ls=linestyles[i], lw=linewidths[i])
		add_trajectory(axs[0], h_rev, c=colors[i], ls=linestyles[i], lw=linewidths[i])
	add_trajectory(axs[0], base_h_fwd,lw=2)
	add_trajectory(axs[0], base_h_rev,lw=2)

	# variation in semi-major axis
	for i in range(4,7):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[0], i, suffix))
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[0], i, suffix))
		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))
			continue
		add_trajectory(axs[1], h_fwd, c=colors[i], ls=linestyles[i], lw=linewidths[i])
		add_trajectory(axs[1], h_rev, c=colors[i], ls=linestyles[i], lw=linewidths[i])
	add_trajectory(axs[1], base_h_fwd,lw=2)
	add_trajectory(axs[1], base_h_rev,lw=2)

	# variation in rotation rate
	for i in range(7,10):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[0], i, suffix))
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[0], i, suffix))
		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))
			continue
		add_trajectory(axs[2], h_fwd, c=colors[i], ls=linestyles[i], lw=linewidths[i])
		add_trajectory(axs[2], h_rev, c=colors[i], ls=linestyles[i], lw=linewidths[i])
	add_trajectory(axs[2], base_h_fwd,lw=2)
	add_trajectory(axs[2], base_h_rev,lw=2)

	# variation in initial time
	for i in range(10,12):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[0], i, suffix))
			add_trajectory(axs[3], h_fwd, c=colors[i], ls=linestyles[i], lw=linewidths[i])
		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))

		try:
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[0], i, suffix))
			add_trajectory(axs[3], h_rev, c=colors[i], ls=linestyles[i], lw=linewidths[i])
		except Exception as ex:
			print("Failed to load {}_orb{}".format(base_star[0], i))
	add_trajectory(axs[3], base_h_fwd,lw=2)
	add_trajectory(axs[3], base_h_rev,lw=2)

	# variation in mass
	star_ind=1
	for i in range(12,14):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[star_ind], i, suffix))
			fwd_diff=(h_fwd.a.values[0]-base_h_fwd.a.values[0])
			add_trajectory(axs[4], h_fwd, c=colors[i], ls=linestyles[i], lw=linewidths[i],a_set=fwd_diff)
		except Exception as ex:
			print(ex)
			print("Failed to load {}_orb{}".format(base_star[star_ind], i))

		try:
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[star_ind], i, suffix))
			rev_diff=(h_rev.a.values[0]-base_h_rev.a.values[0])
			add_trajectory(axs[4], h_rev, c=colors[i], ls=linestyles[i], lw=linewidths[i], a_set=rev_diff)
		except Exception as ex:
			print(ex)
			print("Failed to load {}_orb{}".format(base_star[star_ind], i))

		star_ind+=1
	add_trajectory(axs[4], base_h_fwd,lw=2)
	add_trajectory(axs[4], base_h_rev,lw=2)

	# variation in metallicity
	for i in range(14,16):
		print(i)
		try:
			h_fwd=pd.read_csv("{}/{}_orb{}_{}_fwd/orbital_history.data".format(base_diname, base_star[star_ind], i, suffix))
			add_trajectory(axs[5], h_fwd, c=colors[i], ls=linestyles[i], lw=linewidths[i])
		except Exception as ex:
			print(ex)
			print("Failed to load {}_orb{}".format(base_star[star_ind], i))

		try:
			h_rev=pd.read_csv("{}/{}_orb{}_{}_rev/orbital_history.data".format(base_diname, base_star[star_ind], i, suffix))
			add_trajectory(axs[5], h_rev, c=colors[i], ls=linestyles[i], lw=linewidths[i])
		except Exception as ex:
			print(ex)
			print("Failed to load {}_orb{}".format(base_star[star_ind], i))
			continue

		star_ind+=1
	add_trajectory(axs[5], base_h_fwd,lw=2)
	add_trajectory(axs[5], base_h_rev,lw=2)

	plt.subplots_adjust(hspace=0)

	plt.xlim(-1500,1700)
	plt.show()

def add_trajectory(axs, data, c='k', alpha=1.0, lw=1, ls='-', a_set=0):
	time_normalized=(data.time.values-data.time.values[-1])/1e6
	axs[0].plot(time_normalized, data.e.values, c=c, alpha=alpha, lw=lw, linestyle=ls)
	axs[1].plot(time_normalized, data.a.values-a_set, c=c, alpha=alpha, lw=lw, linestyle=ls)

	axs[0].scatter(time_normalized[-1], data.e.values[-1], c=c, marker='x')
	axs[1].scatter(time_normalized[-1], data.a.values[-1]-a_set, c=c, marker='x')

	axs[0].set_ylabel("e",fontsize=14)
	axs[1].set_ylabel("a [au]",fontsize=14)

	axs[0].set_xticks([-1500,-1000,-500,0,500,1000,1500])
	axs[1].set_xticks([-1500,-1000,-500,0,500,1000,1500])

if __name__=="__main__":
	main()
