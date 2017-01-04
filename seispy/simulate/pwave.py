from __future__ import division
import numpy as np
from ..trace import Trace
from gwpy.frequencyseries import FrequencySeries

def simulate_pwave(stations, A, phi, theta, frequency, duration, Fs=100, c=3000,
        noise_amp=0, phase=0):
    """
    simulate s-wave in a certain direction

    Parameters
    ----------
    direction : TODO
    frequency : TODO
    time : TODO
    Fs : TODO

    Returns
    -------
    E : TODO
    N : TODO
    Z : TODO
    """
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    src_dir = np.array([cphi*stheta, sphi*stheta, ctheta])
    # get time delays
    taus = [np.dot(src_dir, stations[key])/c for key in stations.keys()]
    ts = min(taus)
    te = max(taus)
    times = np.arange(0, ts + duration + te, 1/Fs)
    # shift backward in time
    times -= ts
    data = {}
    for key in stations.keys():
        data[key]={}
        station = stations[key]
        delay = np.dot(src_dir, station)/c
        delaySamps = int(np.round(delay*Fs))
        signal = np.zeros(times.size)
        if frequency == 0:
            signal = A*np.random.randn(times.size)
        else:
            signal = A * np.cos(2*np.pi*frequency*times + phase)
        # impose time delay
        amp = np.roll(signal,delaySamps)
        amp += noise_amp * np.random.randn(times.size)
        data[key]['E'] = Trace(src_dir[0]*amp, sample_rate=Fs, times=times)
        data[key]['N'] = Trace(src_dir[1]*amp, sample_rate=Fs, times=times)
        data[key]['Z'] = Trace(src_dir[2]*amp, sample_rate=Fs, times=times)
        for key2 in data[key].keys():
            data[key][key2].location = station
    return data
