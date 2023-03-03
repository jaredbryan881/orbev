import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

def main():
	fidir="/home/jared/MIT/astero/gyre_HATP2/orbev/data/live_fwd_ps_orot/M1.36_Z0.02262_orb3_live_fwd/"
	finame="orbital_history.data"

	data=pd.read_csv(fidir+finame)
	time=data["time"].values
	a=data["a"].values
	e=data["e"].values
	OmegaRot=data["OmegaRot"].values

	fig,axs=plt.subplots(4,1,sharex=True)
	add_trajectory(time, a, e, OmegaRot, axs)

	axs[3].set_yscale("log")
	axs[3].set_xlabel("Time [Gyr]", fontsize=12)
	axs[0].set_ylabel(r"$a$", fontsize=12)
	axs[1].set_ylabel(r"$e$", fontsize=12)
	axs[3].set_ylabel(r"$dt$ [yr]")
	axs[2].set_ylabel(r"$\Omega_{rot}$", fontsize=12)

	plt.tight_layout()

	plt.show()

def add_trajectory(time, a, e, OmegaRot, axs, color='k'):
	axs[0].plot(time, a, color=color)
	axs[1].plot(time, e, color=color)
	axs[2].plot(time, OmegaRot, color=color)
	axs[3].scatter((time[1:]+time[:-1])/2, time[1:]-time[:-1], color=color, s=2)

if __name__=="__main__":
	main()