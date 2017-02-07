from __future__ import division
import numpy as np
from ..trace import Trace
from gwpy.frequencyseries import FrequencySeries
from astropy import units as u
from numpy.fft.fftpack import ifft, irfft

def simulate_pwave(stations, A, phi, theta, frequency, duration, Fs=100, c=3000,
        noise_amp=0, phase=0, segdur=None):
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
    taus = np.array([-np.dot(src_dir, stations[key])/c for key in
        stations.keys()])
    tau_round = np.round(taus*Fs)/Fs
    ts = min(-tau_round)
    te = max(-tau_round)
    times = np.arange(0, ts + duration + te, 1/Fs)
    Nsamps = duration * Fs
    # shift backward in time
    times += ts
    data = {}
    ct = 0
    final_times = np.arange(0, duration, 1/Fs)
    for key in stations.keys():
        data[key]={}
        station = stations[key]
        delay = -np.dot(src_dir, station)/c
        delaySamps = int(ts*Fs+np.round(delay*Fs))
        signal = np.zeros(times.size)
        if frequency == 0:
            signal = A*np.random.randn(times.size)
        else:
            # most of our noise spectra will be one-sided, but this is a real
            # signal, so we multiply this by two.
            signal = A * np.cos(2*np.pi*frequency*times + phase)
        # impose time delay
        amp = np.roll(signal,delaySamps)[:Nsamps]
        data[key]['E'] = Trace(src_dir[0]*amp, sample_rate=Fs,
                times=final_times, unit=u.m) + gaussian_noise(noise_amp**2, Fs,
                        duration, segdur=segdur)
        data[key]['N'] = Trace(src_dir[1]*amp, sample_rate=Fs,
                times=final_times,unit=u.m) + gaussian_noise(noise_amp**2, Fs,
                        duration, segdur=segdur)
        data[key]['Z'] = Trace(src_dir[2]*amp, sample_rate=Fs,
                times=final_times,unit=u.m) + gaussian_noise(noise_amp**2, Fs,
                        duration, segdur=segdur)
        for key2 in data[key].keys():
            data[key][key2].location = station
    return data

def gaussian_noise(psd_amp, sample_rate, duration, segdur=None):
    if segdur is None:
        segdur=duration
    Nsegs = int(duration / segdur)
    if not(segdur*Nsegs == duration):
        raise('ValueError: segdur must divide duration')
    for seg in range(Nsegs):
        nsamps = sample_rate*segdur
        spec = np.sqrt(2*psd_amp) *\
            np.exp(-1j*2*np.pi*np.random.randn(sample_rate*segdur/2-1))
        if seg ==0:
            ts =\
                np.real(np.fft.ifft(np.concatenate((np.zeros(1),spec,np.zeros(1),np.flipud(spec)*nsamps))))
        else:
            ts = np.hstack((ts,np.real(np.fft.ifft(np.concatenate((np.zeros(1),spec,np.zeros(1),np.flipud(spec)*nsamps))))))

    return Trace(ts, sample_rate=sample_rate, unit=u.m)
