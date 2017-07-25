#! /usr/bin/python
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pymultinest
from utils import *
from seispy.station.stationdata import SeismometerArray
from seispy.station import StationArray, spiral, homestake
from seispy.seispy_io import read_config, print_params
from scipy.linalg import pinv2, svd, inv
from scipy.sparse.linalg import lsqr
import astropy.units as u
import os
import numpy as np
import ConfigParser
import optparse
from collections import OrderedDict
import seispy.plotter.plot as hplot
import matplotlib.tri as mtri
from mpl_toolkits import mplot3d
import healpy as hp
import NX01_AnisCoefficients as anis
import AnisCoefficients_pix as pixAnis

amplitude = 1e-4

# set up array
stations = spiral(50, radius=500)
#stations = homestake('underground')
distances = []
# get some white data
data = SeismometerArray.initialize_all_good(stations, 2000, chans_type='fast')
# get spot size
snames = data.keys()
for station in snames:
    for station2 in snames:
        diff_dir = np.asarray(stations[station]) -\
            np.asarray(stations[station2])
        distances.append(np.sqrt(np.dot(diff_dir,diff_dir)))
maxd = np.max(distances)
print 'Array size is roughly: %d meters' % int(maxd)
spot_size = 5000 / (2 * maxd * 1) / np.sqrt(8)
npix = (4 * np.pi / (spot_size **
        2))#*np.sqrt(len(self.keys()))
# get nside for healpy based on npix (taken up to
# the next power of two)
nside = int(2 ** round(np.log2(round((int(npix)+1) **\
        (0.5)))))
print 'nside is: %d' % nside

lmax = int(np.pi / spot_size)+ 1
print 'lmax roughly is: %d' %  int(lmax+1)

# add some white noise
#data.add_white_noise(1e-24, 100)
# add injection
data.add_p_wave(amplitude, np.radians(0), np.radians(30), 1,2000, c=5000)
sta = data.keys()[0]
asd = data[sta]['HHE'].asd(100)
fftg = fftgram(data[sta]['HHE'], stride=100)
#data.add_s_wave(1e-4, np.radians(-60), np.radians(120), 0, 1, 10, c=3000)
# set up string list
strs = ['s1','s2','p']
strs = ['p']
# get the gamma matrices
#G1, phis1, thetas1, shapes1 = data.get_gamma_matrix_healpy('s', stations, 1, 3000)
# get coherences
Ys = data.get_ffts(1, fftlength=100, overlap=50, window='hann')
# combine gamma matrices
#G = np.vstack((G1, G2)).T
vlist = np.arange(40) * 10 + 4800
vlist = [5000]
loglikes = []
D = Ys.flatten()
D = np.mean(Ys, axis=1)
#nsamps = D.size/(150)
#print D
N_inv = 1e12 * np.diag(np.ones(D.size))
N = 1e-12 * np.diag(np.ones(D.size))
for v in vlist:
    print 'Running velocity: %d' % v
    G2, phis2, thetas2, shapes2 = data.get_response_matrix_healpy('p',
            stations,1,v)
    G = G2.T
    G2 = None
#    G = np.tile(G, (nsamps,1))
    print G.shape
    print D.shape
    # combine angles associated with one axies
    # of gamma matrices. keep track of
    # shapes
    #phis = np.hstack((phis1, phis2))
    phis = phis2
    #thetas = np.hstack((thetas1, thetas2))
    thetas = thetas2
    shapes = []
    #shapes.extend(shapes1)
    shapes.extend(shapes2)
    idx0 = 0

    # average over time
    # do matrix work
    X = np.dot(G.T,D)
#    Gamma_pinv = pinv2(Gamma, rcond=1e-3)
    U, s, Vh = svd(G)
    plt.plot(s)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.savefig('singular_values')
    plt.close()
    G_pinv = pinv2(np.dot(G.T, G), rcond=1e-15)
#    G_pinv = inv(np.dot(G.T, G))
    S = np.real(np.dot(G_pinv, X))

#    S = lsqr(np.real(np.dot(G.T, G)), np.real(np.dot(G.T, Y)))
    print S.shape
    print G.shape
#    cov = np.real(np.dot(np.dot(G.T, N_inv), G))
    G_pinv = None
    cov=None
    model = np.real(np.dot(G,S))
    dat = (D - model)
    ll = -0.5 * np.real(np.dot(np.dot(dat.conj().T, N_inv), dat))
    loglikes.append(ll)

#print loglikes
#loglikes = np.asarray(loglikes[:-1])
#loglikes -= max(loglikes)
#print loglikes
#plt.figure()
#plt.plot(vlist[:-1], loglikes)
#plt.savefig('v_post')
#plt.close()

# make plots
for ii,shape in enumerate(shapes):
    # get shape, get the phis and thetas
    # for this particular map.
    shape = shape[0]
    phis_temp = phis[idx0:idx0+shape]
    thetas_temp = thetas[idx0:idx0+shape]
    # Extract this map.
    S_temp = S[idx0:idx0+shape] / 100
    inj_map = np.zeros(S_temp.shape)
    nside = hp.pixelfunc.npix2nside(shape)
    print nside
    idx = hp.pixelfunc.ang2pix(nside, np.radians(120), np.radians(60))
    inj_map[idx] = 1
    P_inj = pixAnis.clmFromMap(inj_map, lmax-1) * 1e-4
    P_inj[0] = 1e-8
    inj_map = pixAnis.mapFromClm(P_inj, nside)
    ax = plt.subplot(111, projection='mollweide')
    hplot.healpix_heatmap(S_temp, cmap='viridis')
    hplot.outline_text(ax)
    ax.grid(True)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    plt.colorbar()
    ax.set_title(r'$%s$-wave recovery, %d' % (strs[ii],v), fontsize=12)
    plt.savefig('%s_recovery' % (strs[ii]))
    plt.close()
    ax = plt.subplot(111, projection='mollweide')
    hplot.healpix_heatmap(inj_map, cmap='viridis')
    hplot.outline_text(ax)
    ax.grid(True)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    plt.colorbar()
    ax.set_title(r'$%s$-wave injection, %d' % (strs[ii],v), fontsize=12)
    plt.savefig('%s_injection' % (strs[ii]))
    plt.close()

    match = np.dot(inj_map, S_temp) / np.sqrt(np.dot(inj_map,inj_map) *
            np.dot(S_temp, S_temp))
    print 'Match parameter: %f' % match
