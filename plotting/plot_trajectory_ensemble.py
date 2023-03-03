import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def main():
	base_fidir="/home/jared/MIT/astero/gyre_HATP2/orbev/data/live_fwd_ps_orot/"
	fidirs=["M1.36_Z0.02262_orb0_live_fwd",
			"M1.36_Z0.02262_orb1_live_fwd",
			"M1.36_Z0.02262_orb2_live_fwd",
			"M1.36_Z0.02262_orb3_live_fwd",
			"M1.36_Z0.02262_orb4_live_fwd",
			"M1.36_Z0.02262_orb5_live_fwd",
			"M1.36_Z0.02262_orb6_live_fwd",
			"M1.36_Z0.02262_orb7_live_fwd",
			"M1.36_Z0.02262_orb8_live_fwd",
			"M1.36_Z0.02262_orb9_live_fwd",
			"M1.36_Z0.02262_orb10_live_fwd",
			"M1.36_Z0.02262_orb11_live_fwd"]
	"""
	base_fidir="/home/jared/MIT/astero/gyre_HATP2/orbev/data/live_fwd_low_orot/"
	fidirs=["M1.36_Z0.02262_orb0",
			"M1.36_Z0.02262_orb1",
			"M1.36_Z0.02262_orb2",
			"M1.36_Z0.02262_orb3",
			"M1.36_Z0.02262_orb4",
			"M1.36_Z0.02262_orb5",
			"M1.36_Z0.02262_orb6",
			"M1.36_Z0.02262_orb7",
			"M1.36_Z0.02262_orb8",
			"M1.36_Z0.02262_orb9",
			"M1.36_Z0.02262_orb10",
			"M1.36_Z0.02262_orb11"]
	"""

	finame="orbital_history.data"

	fig,axs=plt.subplots(4,1,sharex=True)
	for (i,fidir) in enumerate(fidirs):
		data=pd.read_csv(base_fidir+fidir+"/"+finame)
		time=data["time"].values
		a=data["a"]
		e=data["e"]
		OmegaRot=data["OmegaRot"]

		add_trajectory(time, a, e, OmegaRot, axs)#, color=cm.inferno(i/12))

	axs[3].set_xlabel("Time [Gyr]", fontsize=12)
	axs[0].set_ylabel(r"$a$", fontsize=12)
	axs[1].set_ylabel(r"$e$", fontsize=12)
	axs[3].set_ylabel(r"$dt$ [yr]")
	axs[2].set_ylabel(r"$\Omega_{rot}$", fontsize=12)

	plt.tight_layout()
	plt.show()

def add_trajectory(time, a, e, OmegaRot, axs, color='k'):
	axs[0].plot(time, a, color=color, lw=2)
	axs[1].plot(time, e, color=color, lw=2)
	axs[2].plot(time, OmegaRot, color=color, lw=2)
	axs[3].scatter(((time[1:]+time[:-1])/2), time[1:]-time[:-1], color=color, s=5)

if __name__=="__main__":
	main()