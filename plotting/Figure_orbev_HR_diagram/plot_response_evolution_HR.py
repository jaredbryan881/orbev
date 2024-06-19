import mesa_reader as mr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
import scipy.special as ss
import h5py
import string

def main():
	# MESA grids
	Ms = ["1.20", "1.225", "1.25", "1.275", "1.30", "1.325", "1.35", "1.375", "1.40", "1.425", "1.45", "1.475", "1.50"]
	Zs = ["0.02"]

	plot_mode(Ms, Zs, val_type="amp")

def plot_mode(Ms,Zs, val_type="amp", modes=[29], labels=[r"$g_{29}$"], norms=[[-9,-6]]):
	G = 6.67430e-11 # [m3 km-1 s-2]
	Rsun = 695700000 #[m]
	Msun = 1.989e30 # [kg]

	masses=[1.20, 1.225, 1.25, 1.275, 1.30, 1.325, 1.35, 1.375, 1.40, 1.425, 1.45, 1.475, 1.50]
	fig,axs=plt.subplots(1,len(modes),figsize=(8*len(modes),8), sharex=True, sharey=True)
	for (i,m) in enumerate(Ms):
		color=cm.inferno(i/len(Ms))
		for (j,z) in enumerate(Zs):
			MZ = "M{}_Z{}".format(m,z)
			print(MZ)
			with h5py.File("./conv_and_rad/{}/{}.h5".format(MZ,MZ), "r") as hf:
				star_age=hf["star_age"][:]
				Teff=hf["Teff"][:]
				L=hf["L"][:]
				R=hf["radius"][:]

				xi_r=hf["xi_r"][:]
				lag_L=hf["lag_L"][:]

			xi_r=xi_r[:,:,:,:,0] + 1j*xi_r[:,:,:,:,1]
			lag_L=lag_L[:,:,:,:,0] + 1j*lag_L[:,:,:,:,1]

			inc=np.deg2rad(90)
			omega=np.deg2rad(262)

			mode_amp_all=np.zeros((xi_r.shape[0], 100))
			for t in range(xi_r.shape[0]):
				mode_amp_all[t] = np.log10(np.abs(eval_fourier(xi_r[t], lag_L[t], inc, omega)))

			for (k,mode) in enumerate(modes):
				mode_amp=mode_amp_all[:,mode]
				print("{}: {}".format(mode,np.max(mode_amp)))

				n_pts=n_pts=len(mode_amp)

				points = np.array([np.log10(Teff[:n_pts]), np.log10(L[:n_pts])]).T.reshape(-1, 1, 2)
				segments = np.concatenate([points[:-1], points[1:]], axis=1)

				# Create a continuous norm to map from data points to colors
				norm = plt.Normalize(norms[k][0], norms[k][1])
				lc = LineCollection(segments, cmap='inferno', norm=norm)
				# Set the values used for colormapping
				lc.set_array(mode_amp)
				lc.set_linewidth(15)

				if len(modes)==1:
					cur_ax=axs
					last_ax=axs
				else:
					cur_ax=axs[k]
					last_ax=axs[-1]
				line = cur_ax.add_collection(lc)
				if (i==0) and (k==len(modes)-1):
					cbar=fig.colorbar(line, ax=last_ax)
					cbar.set_label(labels[k], fontsize=14)
				if (i==0):
					cur_ax.text(0.02, 0.95, string.ascii_lowercase[k], transform=cur_ax.transAxes, size=20, weight='bold')
					cur_ax.set_title(labels[k], fontsize=20)

	for k in range(len(modes)):
		if len(modes)==1:
			cur_ax=axs
		else:
			cur_ax=axs[k]
		cur_ax.set_xlim(3.86, 3.76)
		cur_ax.set_ylim(0.2, 0.85)
		cur_ax.set_xlabel(r"$\log_{10}(T_{eff})$ [K]", fontsize=14)
		cur_ax.set_ylabel(r"$\log_{10}(L/L_{\odot})$", fontsize=14)
	plt.tight_layout()
	#plt.savefig("g_amp_across_HR.pdf")
	plt.show()

def eval_fourier (xi_r, lag_L, theta_0, phi_0):
	# Initialize the fourier amplitudes array
	A = np.zeros(xi_r.shape[0], dtype=complex)

	# Loop over l, m and k
	for l in range(2, 3):
		i_l = l-2
		for m in range(-l, l+1):
			i_m = m + 2
			for k in range(1, 100+1):
				i_k = k-1

				Del_R = np.sqrt(4.*np.pi) * xi_r[i_k,i_m,i_l]
				Del_L = np.sqrt(4.*np.pi) * lag_L[i_k,i_m,i_l]
				Del_T = 0.25*(Del_L - 2*Del_R)

				I_0 = 1/2
				I_l = 1/8

				Yml = ss.sph_harm(m, l, phi_0, theta_0)

				Rml = (2 + l)*(1 - l)*I_l/I_0 * Yml
				Tml = 4*I_l/I_0 * Yml

				if k == 0:
					if m == 0:
						kappa = 0.5
					elif m >= 1:
						kappa = 1.
					else:
						kappa = 0.
				else:
					kappa = 1.
				# Add the Fourier contribution
				A[i_k] += 2*kappa*(Del_R*Rml + Del_T*Tml)
	
	# Return data
	return A

if __name__=="__main__":
	main()