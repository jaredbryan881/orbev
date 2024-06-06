import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def main():
	Ms=[1.36, 1.36136, 1.3634, 1.3668, 1.3702, 1.3736, 1.394, 1.428, 1.36, 1.36, 1.36, 1.36, 1.36, 1.36, 1.36]
	Zs=[0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02264262, 0.02267655, 0.0227331, 0.02278965, 0.0228462, 0.0231855, 0.023751]

	with h5py.File("../../data/fixed_orbit/tidal_response_history_allstar.h5", "r") as hf:
		e_dot_freq_allstar=hf["e_dot_freq"][:]
		a_dot_freq_allstar=hf["a_dot_freq"][:]
		e=hf["e"][0]
		a=hf["a"][0]
		Omega_orb=hf["Omega_orb"][:]

	n_times=e_dot_freq_allstar.shape[1]

	e_dot_freq_Ms_allstar=e_dot_freq_allstar[:8]
	e_dot_freq_Zs_allstar=np.concatenate([e_dot_freq_allstar[0][np.newaxis,:], e_dot_freq_allstar[8:]])

	a_dot_freq_Ms_allstar=a_dot_freq_allstar[:8]
	a_dot_freq_Zs_allstar=np.concatenate([a_dot_freq_allstar[0][np.newaxis,:], a_dot_freq_allstar[8:]])

	print("{}-{}".format(e,a))

	e_dot_freq_Ms_diff=np.zeros_like(e_dot_freq_Ms_allstar)
	a_dot_freq_Ms_diff=np.zeros_like(a_dot_freq_allstar)
	e_dot_freq_Zs_diff=np.zeros_like(e_dot_freq_Ms_allstar)
	a_dot_freq_Zs_diff=np.zeros_like(a_dot_freq_allstar)
	for i_star in range(e_dot_freq_Ms_diff.shape[0]):
		e_dot_freq_Ms_diff[i_star]=e_dot_freq_Ms_allstar[i_star]-e_dot_freq_Ms_allstar[0]
		a_dot_freq_Ms_diff[i_star]=a_dot_freq_Ms_allstar[i_star]-a_dot_freq_Ms_allstar[0]

		e_dot_freq_Zs_diff[i_star]=e_dot_freq_Zs_allstar[i_star]-e_dot_freq_Zs_allstar[0]
		a_dot_freq_Zs_diff[i_star]=a_dot_freq_Zs_allstar[i_star]-a_dot_freq_Zs_allstar[0]

	plot_orbev_distributions(e_dot_freq_Ms_allstar, a_dot_freq_Ms_allstar, e_dot_freq_Ms_diff, a_dot_freq_Ms_diff, bins=50)#, foname="orb{}_edot_distribution.pdf".format(i_orb))
	plot_orbev_distributions(e_dot_freq_Zs_allstar, a_dot_freq_Zs_allstar, e_dot_freq_Zs_diff, a_dot_freq_Zs_diff, bins=50)#, foname="orb{}_adot_distribution.pdf".format(i_orb))

	# integrate orbital evolution rate deviation by repeated sampling from the distribution
	e_ensemble_Ms, a_ensemble_Ms = integrate_deviation_distribution(e_dot_freq_Ms_diff, a_dot_freq_Ms_diff)#, foname="orb{}_Ms_horizon.pdf".format(i_orb))
	e_ensemble_Zs, a_ensemble_Zs = integrate_deviation_distribution(e_dot_freq_Zs_diff, a_dot_freq_Zs_diff)#, foname="orb{}_Zs_horizon.pdf".format(i_orb))

	# integrate orbital evolution rates and afterward calcualte the differences
	#e_ensemble_Ms, a_ensemble_Ms = integrate_trajectories(e_dot_freq_Ms_allstar, a_dot_freq_Ms_allstar)#, foname="orb{}_Ms_difference.pdf".format(i_orb))
	#e_ensemble_Zs, a_ensemble_Zs = integrate_trajectories(e_dot_freq_Zs_allstar, a_dot_freq_Zs_allstar)#, foname="orb{}_Zs_difference.pdf".format(i_orb))

def integrate_trajectories(e_dot_freq_allstar, a_dot_freq_allstar, e_sig=0.011, a_sig=0.00065, foname=None):
	n_stars=e_dot_freq_allstar.shape[0]
	n_times=e_dot_freq_allstar.shape[1]

	e_dot=e_dot_freq_allstar.sum(axis=(2,3))
	a_dot=e_dot_freq_allstar.sum(axis=(2,3))
	de=np.cumsum(e_dot, axis=1)
	da=np.cumsum(e_dot, axis=1)

	Delta_e=de[1:,:]-de[0,:]
	Delta_a=da[1:,:]-da[0,:]

	fig,axs=plt.subplots(2,1)
	for i_star in range(n_stars-1):
		ind_fgt_e=first_gt(Delta_e[i_star], e_sig)
		ind_fgt_a=first_gt(Delta_a[i_star], a_sig)

		axs[0].plot(np.arange(len(Delta_e[i_star])), Delta_e[i_star], c=cm.inferno(i_star/8))
		axs[1].plot(np.arange(len(Delta_a[i_star])), Delta_a[i_star], c=cm.inferno(i_star/8))

		axs[0].scatter(ind_fgt_e, Delta_e[i_star][ind_fgt_e], c='k')
		axs[1].scatter(ind_fgt_a, Delta_a[i_star][ind_fgt_a], c='k')
	plt.show()

def integrate_deviation_distribution(e_dot_freq_diff, a_dot_freq_diff, n_ens=1000, foname=None, e_sig=0.011, a_sig=0.00065):
	n_stars=e_dot_freq_diff.shape[0]
	n_times=e_dot_freq_diff.shape[1]

	fig,axs=plt.subplots(1,2,sharex=True,sharey=True,figsize=(10,2.5))
	for i_star in range(n_stars):
		if i_star==0:
			continue
		e_dot_diff=e_dot_freq_diff[i_star].sum(axis=(1,2))
		a_dot_diff=a_dot_freq_diff[i_star].sum(axis=(1,2))

		e_dot_diff=e_dot_diff[~np.isnan(e_dot_diff)]
		a_dot_diff=a_dot_diff[~np.isnan(a_dot_diff)]

		signs_e_dot=np.sign(e_dot_diff)
		signs_a_dot=np.sign(e_dot_diff)

		e_dot_diff=np.abs(e_dot_diff)
		a_dot_diff=np.abs(a_dot_diff)

		ntimes_cur=len(e_dot_diff)

		# interpolate differences to a finer grid for integration
		times=np.linspace(2.6e9, 2.6e9+1e5*ntimes_cur, n_times)[:ntimes_cur]

		de_diff=e_dot_diff*1e5
		da_diff=a_dot_diff*1e5

		inds_fgt_e=np.zeros(n_ens)
		inds_fgt_a=np.zeros(n_ens)
		for i_rzn in range(n_ens):
			# mimic "drawing from the ensemble" as simple a shuffling of the discretized edots
			cur_inds=np.random.permutation(ntimes_cur)
			e_ensemble=np.cumsum(de_diff[cur_inds])
			a_ensemble=np.cumsum(da_diff[cur_inds])

			# find spot where e goes above e_sig
			ind_fgt_e=first_gt(e_ensemble, e_sig)
			inds_fgt_e[i_rzn]=ind_fgt_e

			ind_fgt_a=first_gt(a_ensemble, a_sig)
			inds_fgt_a[i_rzn]=ind_fgt_a

		inds_fgt_e=inds_fgt_e[np.nonzero(inds_fgt_e)]*1e5
		inds_fgt_a=inds_fgt_a[np.nonzero(inds_fgt_a)]*1e5

		axs[0].hist(np.log10(inds_fgt_e), bins=20, fc=cm.inferno(i_star/8), alpha=0.25, density=True)
		axs[1].hist(np.log10(inds_fgt_a), bins=20, fc=cm.inferno(i_star/8), alpha=0.25, density=True)

	axs[1].set_xlabel("Time [yr]")
	axs[0].set_ylabel(r"$T_e$ density")
	axs[1].set_ylabel(r"$T_a$ density")

	if foname is not None:
		plt.savefig(foname)
		plt.close()
	else:
		plt.show()

	return e_ensemble, a_ensemble

def first_gt(arr, val):
	return np.argmax(arr>val)

def plot_orbev_distributions(e_dot_freq_allstar, a_dot_freq_allstar, e_dot_freq_diff, a_dot_freq_diff, bins=50, foname=None):
	fig,axs=plt.subplots(2,2, figsize=(10,5), sharey='row', sharex='row')
	for i_star in range(e_dot_freq_allstar.shape[0]):
		e_dot_freq_allstar_sum=e_dot_freq_allstar[i_star].sum(axis=(1,2))
		e_dot_freq_allstar_sum=e_dot_freq_allstar_sum[~np.isnan(e_dot_freq_allstar_sum)]
		e_dot_freq_allstar_sum=np.log10(np.abs(e_dot_freq_allstar_sum))

		width=(e_dot_freq_allstar_sum.max()-e_dot_freq_allstar_sum.min())/bins
		hist, bin_edges = np.histogram(e_dot_freq_allstar_sum, bins=bins, density=True)
		axs[0,0].bar(bin_edges[:-1], hist, color=cm.inferno(i_star/8), alpha=0.25, width=width)

		a_dot_freq_allstar_sum=a_dot_freq_allstar[i_star].sum(axis=(1,2))
		a_dot_freq_allstar_sum=a_dot_freq_allstar_sum[~np.isnan(a_dot_freq_allstar_sum)]
		a_dot_freq_allstar_sum=np.log10(np.abs(a_dot_freq_allstar_sum))

		width=(a_dot_freq_allstar_sum.max()-a_dot_freq_allstar_sum.min())/bins
		hist, bin_edges = np.histogram(a_dot_freq_allstar_sum, bins=bins, density=True)
		axs[0,1].bar(bin_edges[:-1], hist, color=cm.inferno(i_star/8), alpha=0.25, width=width)

		if i_star>0:
			e_dot_freq_diff_sum=e_dot_freq_diff[i_star].sum(axis=(1,2))
			e_dot_freq_diff_sum=e_dot_freq_diff_sum[~np.isnan(e_dot_freq_diff_sum)]
			e_dot_freq_diff_sum=np.log10(np.abs(e_dot_freq_diff_sum))

			width=(e_dot_freq_diff_sum.max()-e_dot_freq_diff_sum.min())/bins
			hist, bin_edges = np.histogram(e_dot_freq_diff_sum, bins=bins, density=True)
			axs[1,0].bar(bin_edges[:-1], hist, color=cm.inferno(i_star/8), alpha=0.25, width=width)

			a_dot_freq_diff_sum=a_dot_freq_diff[i_star].sum(axis=(1,2))
			a_dot_freq_diff_sum=a_dot_freq_diff_sum[~np.isnan(a_dot_freq_diff_sum)]
			a_dot_freq_diff_sum=np.log10(np.abs(a_dot_freq_diff_sum))

			width=(a_dot_freq_diff_sum.max()-a_dot_freq_diff_sum.min())/bins
			hist, bin_edges = np.histogram(a_dot_freq_diff_sum, bins=bins, density=True)
			axs[1,1].bar(bin_edges[:-1], hist, color=cm.inferno(i_star/8), alpha=0.25, width=width)

		axs[0,0].set_yscale("log")
		axs[0,1].set_yscale("log")
		axs[1,0].set_yscale("log")
		axs[1,1].set_yscale("log")
	
	axs[0,0].set_xlabel(r"$\dot{e}$ [1/yr]")
	axs[0,1].set_xlabel(r"$\dot{a}$ [au/yr]")

	axs[1,0].set_xlabel(r"$\Delta \dot{e}$ [1/yr]")
	axs[1,1].set_xlabel(r"$\Delta \dot{a}$ [au/yr]")

	axs[0,0].set_ylabel("Density")
	axs[1,0].set_ylabel("Density")

	plt.tight_layout()

	if foname is not None:
		plt.savefig(foname)
		plt.close()
	else:
		plt.show()

if __name__=="__main__":
	main() 