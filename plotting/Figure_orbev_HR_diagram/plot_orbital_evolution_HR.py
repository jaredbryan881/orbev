import mesa_reader as mr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
import h5py

def main():
	# MESA grids
	Ms = ["1.20", "1.225", "1.25", "1.275", "1.30", "1.325", "1.35", "1.375", "1.40", "1.425", "1.45", "1.475", "1.50"]
	Zs = ["0.02"]

	plot_orbital_evolution(Ms, Zs, val_str="e_dot", label=r"$\log_{10}(|\dot{e}|)$ [1/yr]", clim=[-10.0, -3.0], foname=None)

def plot_orbital_evolution(Ms, Zs, val_str="e_dot", label=r"$\log_{10}(|\dot{e}|)$ [1/yr]", clim=None, foname=None):
	G = 6.67430e-11  # [m3 km-1 s-2]
	Rsun = 695700000 # [m]
	Msun = 1.989e30  # [kg]
	au = 1.496e+11   # [m]

	masses=[1.20, 1.225, 1.25, 1.275, 1.30, 1.325, 1.35, 1.375, 1.40, 1.425, 1.45, 1.475, 1.50]
	fig,axs=plt.subplots(1,1,figsize=(8.5,7))
	for (i,m) in enumerate(Ms):
		color=cm.inferno(i/len(Ms))
		for (j,z) in enumerate(Zs):
			MZ = "M{}_Z{}".format(m,z)
			with h5py.File("./conv_and_rad/{}/{}.h5".format(MZ,MZ), "r") as hf:
				star_age=hf["star_age"][:]
				val=hf[val_str][:]
				Teff=hf["Teff"][:]
				L=hf["L"][:]
				R=hf["radius"][:]

			n_pts=len(val)

			M = masses[i]
			R = R[:n_pts]
			freq_scale = (86400/(2*np.pi))*np.sqrt((G*M*Msun)/((R*Rsun)**3)) # [1/day]

			val=val[:,:,:,0,0]
			val=val.sum(axis=(1,2))

			# dimensionalize the value
			val*=freq_scale*365 # [1/yr]
			# more to be done if it's a_dot, J_dot, or o_dot
			if val_str=="a_dot":
				val*=R*(Rsun/au) # now in au/yr
				val_norm=np.log10(np.abs(val/freq_scale**(5/3)))
				print(np.max(val_norm))
			elif val_str=="e_dot":
				val_norm=np.log10(np.abs(val/freq_scale**(7/3))) # keep in 1/yr
				print(np.max(val_norm))
			elif val_str=="o_dot":
				val_norm=np.rad2deg(val) # now in deg/yr
			elif val_str=="J_dot":
				# TODO: dimensionalize Jdot
				val_norm=np.log10(np.abs(val))
				pass
			else:
				print("how did you even get here")

			points = np.array([np.log10(Teff[:n_pts]), np.log10(L[:n_pts])]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)

			colors=[]
			for v in val:
				pct=(np.log10(np.abs(v))-clim[0])/(clim[1]-clim[0])
				if v>0:
					rgba=cm.Reds(pct)
				elif v<=0:
					rgba=cm.Blues(pct)
				colors.append(rgba)

			# Create a continuous norm to map from data points to colors
			if clim is not None:
				norm = plt.Normalize(clim[0], clim[1])
			else:
				norm = plt.Normalize(val_norm.min(), val_norm.max())
			lc = LineCollection(segments, color=colors, norm=norm)
			# Set the values used for colormapping
			lc = LineCollection(segments, cmap='inferno', norm=norm)
			lc.set_array(val_norm)
			lc.set_linewidth(20)
			line = axs.add_collection(lc)
			if i==0:
				cbar=fig.colorbar(line, ax=axs)
				cbar.set_label(label, fontsize=14)
	axs.set_xlim(3.86, 3.76)
	axs.set_ylim(0.2, 0.85)
	axs.set_xlabel(r"$\log_{10}(T_{eff})$ [K]", fontsize=14)
	axs.set_ylabel(r"$\log_{10}(L/L_{\odot})$", fontsize=14)
	plt.tight_layout()
	if foname is not None:
		plt.savefig("edot_across_HR.pdf", rasterized=True)
		plt.close()
	else:
		plt.show()

if __name__=="__main__":
	main()