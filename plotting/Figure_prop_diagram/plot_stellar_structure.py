import numpy as np
import mesa_reader as mr
import matplotlib.pyplot as plt
import matplotlib.cm as cm

Rsun = 695700000 #[m]

def main():
	base_fidir="/home/jared/MIT/astero/gyre_HATP2/orbev/data/stellar_structure/"

	plot_propagation_diagram(base_fidir)


def plot_propagation_diagram(base_fidir):
	data_ZAMS=mr.MesaData(base_fidir+"ZAMS/profile6.data")
	data_MID=mr.MesaData(base_fidir+"MID/profile6.data")
	data_REMS=mr.MesaData(base_fidir+"REMS/profile6.data")

	fs = (86400/2*np.pi)/1e6

	fig,axs=plt.subplots(3,1, sharex=True, sharey=True, figsize=(8,8))
	# plot ZAMS
	add_propagation_diagram(axs[0], 10**(data_ZAMS.logR), data_ZAMS.brunt_frequency, data_ZAMS.lamb_Sl2*fs)

	# plot MID
	add_propagation_diagram(axs[1], 10**(data_MID.logR), data_MID.brunt_frequency, data_MID.lamb_Sl2*fs)

	# plot REMS
	add_propagation_diagram(axs[2], 10**(data_REMS.logR), data_REMS.brunt_frequency, data_REMS.lamb_Sl2*fs)
	

	axs[2].set_xlabel(r"$Radius$ $[R_{\odot}]$", fontsize=12)
	axs[0].set_ylabel(r"$log_{10}(\sigma$ [cyc/day] $)$", fontsize=12)
	axs[1].set_ylabel(r"$log_{10}(\sigma$ [cyc/day] $)$", fontsize=12)
	axs[2].set_ylabel(r"$log_{10}(\sigma$ [cyc/day] $)$", fontsize=12)

	axs[0].set_yscale("log")
	axs[0].set_ylim(1e-1, 1e3)
	axs[0].set_xlim(0,1.95)

	plt.tight_layout()
	plt.show()

def add_propagation_diagram(ax, radius, N, Sl2):
	# Brunt-Vaisala Frequency [cyc/day]
	ax.plot(radius, N, c='steelblue', lw=2)

	# l=2 Lamb frequency [cyc/day]
	ax.plot(radius, Sl2, c='crimson', lw=2)

	# g-modes
	ax.fill_between(radius, N, np.ones(len(radius))*1e-3, alpha=0.5, color='steelblue')

	# p-modes
	ax.fill_between(radius, Sl2, np.ones(len(radius))*1e5, alpha=0.5, color='crimson')

if __name__=="__main__":
	main()