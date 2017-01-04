from __future__ import division
import numpy as np
from seispy.utils import orf_p_directional
from seispy.simulate import simulate_pwave
from scipy.sparse.linalg import lsqr

stations = {}

stations['1']=[0,0,0]
stations['2']=[0, 0, 1000]
stations['3']=[1000, 0, 0]
stations['4']=[0,1000,0]


data = simulate_pwave(stations, 1, 0, np.pi/2, 1, 10000)
First = 1
for ii in '1234':
    for jj in '1234':
        if int(ii)>=int(jj):
            continue
        else:
            p1 = data[ii]['Z'].average_fft(fftlength=10)
            p2 = data[jj]['Z'].average_fft(fftlength=10)
            idx = np.where(p1.frequencies.value==1)
            p12 = np.conj(p1[idx])*p2[idx]
            p1 = 1
            p2 = 1
            Y.append(spec[idx])
        gammas, phis, thetas = orf_p_directional([0,0,1],[0,0,1],stations['1'], stations['2'], 3000, 1,
            thetamesh=6, phimesh=6)
        if First == 1:
            GG = np.dot(np.conj(gamma), np.transpose(gamma))
            GY = np.conj(gamma)*p12
            First = 0
        else:
            GG += np.dot(np.conj(gamma), np.transpose(gamma))
            GY = np.conj(gamma)*p12
S = lsqr(np.real(GG), np.real(GY), 1e-6, 1000)
plt.figure()
plt.subplot(111,projection='aitoff')
plt.pcolormesh(phis-np.pi, thetas-np.pi/2, S, 'viridis')
