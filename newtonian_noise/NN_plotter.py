from matplotlib import pyplot as plt
import healpy as hp
import numpy as np

def plot_ampl_map(coherent_map,pol,filename=None):
	thetas=coherent_map['thetas']
	phis=coherent_map['phis']
	ampl=coherent_map.get_amplitude(pol)

	hp.mollview(ampl)
	hp.graticule()

	if filename is None:
		plt.show()
	else:
		plt.savefig(filename)