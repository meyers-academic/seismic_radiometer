from __future__ import division
import numpy as np
from ..trace import Trace
from gwpy.frequencyseries import FrequencySeries

def simulate_swave(stations, A, phi, theta, psi, frequency, duration, Fs=100, c=3000,
        noise_amp=0, phase=0):
    """
    simulate s-wave in a certain direction

    Parameters
    ----------
    stations : `dict`
        dictionary of station locations
    A : `float`
        amplitude of input wave
    phi : `float`
        azimuth in radians
    theta : `float`
        polar angle from north pole in radians
    psi : `float`
        s-wave polarization angle from horizontal
        E-N plane in radians
    frequency : `float`
        frequency of source
    duration : `float`
        duration of signal to simulate
    Fs : `float`, optional, default=100 Hz
        sample rate (int preferred)
    c : `float`, optional, default=3000 m/s
        speed of wave
    noise_amp : `float`, optional, default=0
        noise amplitude (in m^2/Hz)
    phase : `float`, optional, default=0
        phase delay of wave in radians

    Returns
    -------
    data : `dict`
        2-layer dict with first keys as stations,
        second keys as channels for each station.
        Each entry is the data for that channel
        for that station for a simulated wave.
    """
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    src_dir = np.array([cphi*stheta, sphi*stheta, ctheta])
    # Get relative amplitudes in E,N,Z directions
    # based on polarizations. See internal method below.
    dx, dy, dz = get_polarization_coeffs(phi, theta, psi)

    # get time delays
    taus = np.array([-np.dot(src_dir, stations[key])/c for key in
        stations.keys()])
    tau_round = np.round(taus*Fs)/Fs
    ts = min(-tau_round)
    te = max(-tau_round)
    Nsamps = duration * Fs
    final_times = np.arange(0, duration, 1/Fs)
    times = np.arange(0, ts + duration + te, 1/Fs)
    # shift backward in time
    times += ts
    data = {}
    ct = 0
    for key in stations.keys():
        data[key]={}
        station = stations[key]
        delay = -np.dot(src_dir, station)/c
        delaySamps = int(ts*Fs+np.round(delay*Fs))
        signal = np.zeros(times.size)
        if frequency == 0:
            signal = A*np.random.randn(times.size)
        else:
            signal = A * np.sin(2*np.pi*frequency*times + phase)
        # impose time delay
        amp = np.roll(signal,delaySamps)[:Nsamps]
        amp += noise_amp * np.random.randn(Nsamps)
        data[key]['E'] = Trace(dx*amp, sample_rate=Fs, times=final_times)
        data[key]['N'] = Trace(dy*amp, sample_rate=Fs, times=final_times)
        data[key]['Z'] = Trace(dz*amp, sample_rate=Fs, times=final_times)
        for key2 in data[key].keys():
            data[key][key2].location = station
    return data

def get_polarization_coeffs(phi, theta, psi):
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    # Setting up polarization coefficients is tricky.
    # Decide horizontal is psi=0.
    # polarization bases are:
    # p1 = [-cphi*ctheta, -sphi*ctheta, stheta]
    # p2 = [-sphi, cphi, 0]
    # pol_coeffs = cpsi*p2 + spsi*p1
    # dx, dy, dz will be pol_coeffs.
    dx = -ctheta*cphi*spsi - cpsi*sphi
    dy = -spsi*ctheta*sphi + cpsi*cphi
    dz = spsi*stheta
    return dx, dy, dz

