#! /usr/bin/python
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pymultinest
from seispy.station.stationdata import SeismometerArray
from seispy.station import StationArray, spiral, homestake
from seispy.seispy_io import read_config, print_params
from scipy.linalg import pinv2,svd
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
stations = spiral(9, radius=500)
#stations = homestake('underground')
distances = []
# get some white data
data = SeismometerArray.initialize_all_good(stations, 1000, chans_type='fast')
# get spot size
snames = data.keys()
for station in snames:
    for station2 in snames:
        diff_dir = np.asarray(stations[station]) -\
            np.asarray(stations[station2])
        distances.append(np.sqrt(np.dot(diff_dir,diff_dir)))

spot_size = 5000 / (2 * 1 * max(distances))
lmax = int(np.pi / spot_size)+ 1
print 'lmax roughly is: %d' %  int(lmax+1)

# add some white noise
#data.add_white_noise(1e-24, 100)
# add injection
data.add_p_wave(amplitude, np.radians(-120), np.radians(-120), 1,1000, c=5000)
data.add_s_wave(1e-6, np.radians(-60), np.radians(120), 0, 1, 1000, c=3000)
# set up string list
strs = ['s1','s2','p']
#strs = ['p']
# get the gamma matrices
#G1, phis1, thetas1, shapes1 = data.get_gamma_matrix_healpy('s', stations, 1, 3000)
# get coherences
Ys = data.get_coherences(1, fftlength=100, overlap=50, window='hann')
# combine gamma matrices
#G = np.vstack((G1, G2)).T
vlist = np.arange(40) * 10 + 4800
vlist = [5000]#, 5000]
loglikes = []
Y = np.mean(Ys, 1)
print Y.size
N_inv = 1e24 * np.diag(np.ones(Y.size))
N = 1e-24 * np.diag(np.ones(Y.size))
nsamps = 10
for v in vlist:
    print 'Running velocity: %d' % v
    G2, phis2, thetas2, shapes2 = data.get_gamma_matrix_healpy('p',
            stations,1,v)
    G1, phis1, thetas1, shapes1 = data.get_gamma_matrix_healpy('s', stations,
            1, 3000)
    shapes = []
    shapes.extend(shapes1)
    shapes.extend(shapes2)
    G = np.vstack((G2, G1)).T
    G2 = None
    G1 = None
    #G = np.tile(G, (nsamps,1))
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
    # get SVD
    U,s,Vh = svd(G)
    rank = np.sum(s > 1e-5 * np.max(s))
    print 'Rank of matrix is: %d' % rank
    plt.figure()
    plt.plot(s)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.savefig('svd_test')
    # average over time
    # do matrix work
    X = np.real(np.dot(G.T,Y))
#    Gamma_pinv = pinv2(Gamma, rcond=1e-3)
    G_pinv = pinv2(np.real(np.dot(G.T, G)), rcond=1e-5)
    S = np.real(np.dot(G_pinv, X))

#    S = lsqr(np.real(np.dot(G.T, G)), np.real(np.dot(G.T, Y)))
#    cov = np.real(np.dot(np.dot(G.T, N_inv), G))
    cov=None
    model = np.real(np.dot(G,S))
    dat = (Y - model)
    ll = -0.5 * np.real(np.dot(np.dot(dat.conj().T, N_inv), dat))
    loglikes.append(ll)
    # Model resolution matrix
    MRM = np.dot(G_pinv, np.real(np.dot(G.T, G)))
    G_pinv = None
    plt.figure()
    plt.pcolormesh(MRM)
    plt.colorbar()
    plt.savefig('model_resolution_matrix')
    plt.close()
    MRM_diag = np.diag(MRM)

# make plots
for ii,shape in enumerate(shapes):
    # get shape, get the phis and thetas
    # for this particular map.
    print shape
    shape = shape[0]
    phis_temp = phis[idx0:idx0+shape]
    thetas_temp = thetas[idx0:idx0+shape]
    # Extract this map.
    S_temp = S[idx0:idx0+shape] / 100
    inj_map = np.zeros(S_temp.shape)
    nside = 8
    idx = hp.pixelfunc.ang2pix(nside, np.radians(120), np.radians(60))
    inj_map[idx] = 1
    P_inj = pixAnis.clmFromMap(inj_map, lmax+1) * 1e-8
    P_inj[0] = 1e-12
    inj_map = pixAnis.mapFromClm(P_inj, nside)
    # plot recovery map
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
    # plot injection map
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
    # plot model resolution matrix diagonals
    ax = plt.subplot(111, projection='mollweide')
    hplot.healpix_heatmap(MRM_diag, cmap='viridis')
    hplot.outline_text(ax)
    ax.grid(True)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    plt.colorbar()
    ax.set_title(r'$%s$-wave model resolution, %d' % (strs[ii],v), fontsize=12)
    plt.savefig('%s_model_resolution' % (strs[ii]))
    plt.close()




    match = np.dot(inj_map, S_temp) / np.sqrt(np.dot(inj_map,inj_map) *
            np.dot(S_temp, S_temp))
    print 'Match parameter: %f' % match
