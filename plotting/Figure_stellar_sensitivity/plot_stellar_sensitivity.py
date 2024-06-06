import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def main():
	Ms=[1.36, 1.36136, 1.3634, 1.3668, 1.3702, 1.3736, 1.394, 1.428, 1.36, 1.36, 1.36, 1.36, 1.36, 1.36, 1.36]
	Zs=[0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02262, 0.02264262, 0.02267655, 0.0227331, 0.02278965, 0.0228462, 0.0231855, 0.023751]

	dt=1e5

	fig1,axs1=plt.subplots(8,2,sharex=True, sharey=True, figsize=(10,10))
	fig2,axs2=plt.subplots(8,2,sharex=True, sharey=True, figsize=(10,10))
	for (i, M) in enumerate(Ms):
		try:
			with h5py.File("../../data/fixed_orbit/tidal_response_history_{}.h5".format(i)) as hf:
				e_dot_freq=hf["e_dot_freq"][:]
		except:
			continue

		M_cur=Ms[i]
		Z_cur=Zs[i]
		fs=np.loadtxt("../../data/fixed_orbit/fs_M{}_Z{}.txt".format(M_cur, Z_cur))
		fs=fs[:e_dot_freq.shape[0]]

		npts=e_dot_freq.shape[0]
		e_dot_total=np.sum(e_dot_freq, axis=2)
		e_dot_total[:,0]*=fs
		e_dot_total[:,1]*=fs 
		e_dot_total*=365 # 1/yr

		if i==0:
			# Mass perturbations
			axs1[i,0].scatter((2.6e9 + np.arange(npts)*dt)/1e6, np.abs(e_dot_total[:,0]), c='k', s=1)
			axs1[i,0].set_yscale("log")

			axs1[i,1].scatter((2.6e9 + np.arange(npts)*dt)/1e6, np.abs(e_dot_total[:,1]), c='k', s=1)
			axs1[i,1].set_yscale("log")

			# Metallicity perturbations
			axs2[i,0].scatter((2.6e9 + np.arange(npts)*dt)/1e6, np.abs(e_dot_total[:,0]), c='k', s=1)
			axs2[i,0].set_yscale("log")

			axs2[i,1].scatter((2.6e9 + np.arange(npts)*dt)/1e6, np.abs(e_dot_total[:,1]), c='k', s=1)
			axs2[i,1].set_yscale("log")

			axs1[i,0].set_ylabel(r"$\dot{e}$ [1/yr]", fontsize=12)
			axs2[i,0].set_ylabel(r"$\dot{e}$ [1/yr]", fontsize=12)

		elif i<8:
			# Mass perturbations
			axs1[i,0].scatter((2.6e9 + np.arange(npts)*dt)/1e6, np.abs(e_dot_total[:,0]), c='k', s=1)
			axs1[i,0].set_yscale("log")

			axs1[i,1].scatter((2.6e9 + np.arange(npts)*dt)/1e6, np.abs(e_dot_total[:,1]), c='k', s=1)
			axs1[i,1].set_yscale("log")

			axs1[i,0].set_title("M+{:.2f}%".format(100*(Ms[i]-Ms[0])/Ms[0]))
			axs1[i,1].set_title("M+{:.2f}%".format(100*(Ms[i]-Ms[0])/Ms[0]))

			axs1[i,0].set_ylabel(r"$\dot{e}$ [1/yr]", fontsize=12)

		elif i>=8:
			# Metallicity perturbations
			axs2[i-7,0].scatter((2.6e9 + np.arange(npts)*dt)/1e6, np.abs(e_dot_total[:,0]), c='k', s=1)
			axs2[i-7,0].set_yscale("log")

			axs2[i-7,1].scatter((2.6e9 + np.arange(npts)*dt)/1e6, np.abs(e_dot_total[:,1]), c='k', s=1)
			axs2[i-7,1].set_yscale("log")

			axs2[i-7,0].set_title("Z+{:.2f}%".format(100*(Zs[i]-Zs[0])/Zs[0]))
			axs2[i-7,1].set_title("Z+{:.2f}%".format(100*(Zs[i]-Zs[0])/Zs[0]))

			axs2[i-7,0].set_ylabel(r"$\dot{e}$ [1/yr]", fontsize=12)

	axs1[-1,0].set_xlabel("Time [Myr]", fontsize=12)
	axs1[-1,1].set_xlabel("Time [Myr]", fontsize=12)

	axs1[0,0].set_title("m=0",fontsize=12)
	axs1[0,1].set_title("m=2",fontsize=12)

	axs2[-1,0].set_xlabel("Time [Myr]", fontsize=12)
	axs2[-1,1].set_xlabel("Time [Myr]", fontsize=12)

	axs2[0,0].set_title("m=0",fontsize=12)
	axs2[0,1].set_title("m=2",fontsize=12)
	
	fig1.tight_layout()
	fig2.tight_layout()

	#fig1.savefig("edot_Ms_orb{}.png".format(i_orb))
	#fig2.savefig("edot_Zs_orb{}.png".format(i_orb))

	#plt.close()
	plt.show()

if __name__=="__main__":
	main()