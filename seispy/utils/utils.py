from __future__ import division
import numpy as np
from matplotlib import use,rc
use('agg')
rc('text',usetex=True)
from seispy.trace import Trace
import astropy.units as u
import os
from lal import CreateREAL8FrequencySeries
from lal import CreateREAL8TimeSeries
from lalsimulation import SimNoise
from gwpy.spectrogram import Spectrogram
import pickle

# save object
def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

# load object
def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def fftgram(ts, stride):
    """Calculate the Fourier-gram of this `TimeSeries`.
    At every ``stride``, a single, complex FFT is calculated.
    Parameters
    ----------
    stride : `float`
        number of seconds in single PSD (column of spectrogram)
    Returns
    -------
    fftgram : :class:`~gwpy.spectrogram.core.Spectrogram`
        a Fourier-gram
    """
    fftlength = stride
    dt = stride
    df = 1/fftlength
    stride *= ts.sample_rate.value
    # get size of Spectrogram
    nsteps = int(ts.size // stride)
    # get number of frequencies
    nfreqs = int(fftlength*ts.sample_rate.value/2)

    # generate output spectrogram
    dtype = np.complex
    dat = np.zeros((nsteps, nfreqs), dtype=dtype)
    # stride through TimeSeries, recording FFTs as columns of Spectrogram
    for step in range(nsteps):
        # find step TimeSeries
        idx = int(stride * step)
        idx_end = int(idx + stride)
        stepseries = ts[idx:idx_end]
        # calculated FFT and stack
        stepfft = stepseries.fft()
        dat[step] = stepfft.value[1:]
    out = Spectrogram(dat,
                      name=ts.name, t0=ts.t0, f0=0, df=df,
                      dt=dt, copy=False, unit=ts.unit*(u.Hz**-0.5), dtype=dtype)

    return out

def set_channel_vector(letter):
    if letter.upper() == 'Z':
        return [0, 0, 1]
    if letter.upper() == 'N':
        return [0, 1, 0]
    if letter.upper() == 'E':
        return [1, 0 ,0]
    if letter.upper() == 'HHZ':
        return [0, 0, 1]
    if letter.upper() == 'HHN':
        return [0, 1, 0]
    if letter.upper() == 'HHE':
        return [1, 0 ,0]

def _calc_travel_time(delta_vec, Omega, v):
    """
    calculate travel time between two stations assuming
    the ray takes a straight line.
    """
    delta_vec = np.array(delta_vec)
    Omega = np.array(Omega)
    dt = np.dot(delta_vec, Omega)/v
    return dt

def calc_travel_time(delta_vec, OMEGA, v):
    """
    Returns light travel time for many directions on the sky.

    Parameters
    ----------
    delta_vec : `numpy.ndarray`, list
        vector for detector separation. n x 1 matrix.
    OMEGA : `numpy.ndarray`
        matrix of sky direction unit vectors. m x n matrix.
        m is number of sky directions, n is dimension of space.
    v : `float`
        velocity of seismic wave

    Returns
    -------
    dt : `numpy.ndarray`
        list of travel times for each sky direction
    """
    delta_vec = np.array(delta_vec)
    OMEGA = np.array(OMEGA)
    omg_shape = OMEGA.shape
    if not(delta_vec.size == omg_shape[1]):
        raise ValueError('Delta vec should have same dimension as columns of\
        Omega vec')
    dt = np.array([_calc_travel_time(delta_vec, OMEGA[ii],v) for ii in
        range(omg_shape[0])])
    return dt

def calc_travel_time2(delta_vec, OMEGA, v):
    delta_vec = np.array(delta_vec)
    OMEGA = np.array(OMEGA)
    return np.dot(delta_vec, OMEGA.T)/v

def get_pdf_from_map(M):
    """
    get pdf for joint phi/theta and each
    individually from a map

    Parameters
    ----------
    M : TODO

    Returns
    -------
    TODO

    """
    M[M<0] = 0
    pdf = M/M.flatten().sum()
    # marg over theta
    pdf_phi = pdf.sum(axis=1)
    # marg over phi
    pdf_theta = pdf.sum(axis=0)
    return pdf, pdf_phi, pdf_theta

def get_1d_conf(pdf, vals, conf=0.90):
    """
    get 1d confidence interval using descending probabilities method

    Parameters
    ----------
    pdf : `numpy.ndarray`
        posterior to get confidence interval for

    Returns
    -------
    ci : `tuple`
        (low, high) ends of confidence interval
    """
    idxs = np.argsort(np.asarray(pdf))[::-1]
    cdf = pdf[idxs].cumsum()
    conf_args = idxs[cdf < conf]
    low = min(vals[conf_args])
    high = max(vals[conf_args])
    return (low,high)

def get_phi_theta_intervals_from_map(M, phis, thetas, conf=0.90):
    """
    get phi and theta confidence intervals from map

    Parameters
    ----------
    M : `numpy.ndarray`
        recovery map
    conf : `float`
        confidence interval value

    Returns
    -------
    ci_phi : `tuple`
        phi confidence interval
    ci_theta : `tumple`
        theta confidence interval
    """
    pdf, pdf_phi, pdf_theta = get_pdf_from_map(M)
    ci_phi = get_1d_conf(pdf_phi, phis, conf=conf)
    ci_theta = get_1d_conf(pdf_theta, thetas, conf=conf)
    return ci_phi, ci_theta

def combine_data_dicts(dict_list):
    """
    combine a list of data dicts into one final list.

    Parameters
    ----------
    dict_list : `list` of `dicts`
        list of data dictionaries

    Returns
    -------
    final_dict : `dict`
        final data dictionary
    """
    First = 1
    for D in dict_list:
        if First:
            final_dict = D
            First = 0
        else:
            final_dict = _combine_data_dicts(D, final_dict)
    return final_dict

def _combine_data_dicts(dict1,dict2):
    """combine two data dicts together and return a third

    Parameters
    ----------
    dict1 : TODO
    dict2 : TODO

    Returns
    -------
    TODO

    """
    dict3 = {}
    for station in dict1.keys():
        dict3[station] = {}
        for channel in dict1[station].keys():
            tr1 = dict1[station][channel]
            tr2 = dict2[station][channel]
            Fs = tr1.sample_rate.value
            # check which stream starts first
            d1_first = tr1.times.value[0] > tr2.times.value[0]
            # if it's trace 1
            start_diff = tr1.times.value[0] - tr2.times.value[0]
            end_diff = tr1.times.value[-1] - tr2.times.value[-1]
            start_zeros = np.zeros(np.abs(start_diff)*tr1.sample_rate.value)
            if np.abs(end_diff) > 0:
                end_zeros = np.zeros(np.abs(end_diff)*tr1.sample_rate.value - 1)
            # trace 1 starts after trace 2. Add zeros to beginning of trace 1
            # if trace 1 ends after trace 2, add zeros to end of trace 2
            if start_diff==0 and end_diff==0:
                tr1_new = tr1.copy()
                tr2_new = tr2.copy()
            elif start_diff >= 0 and end_diff >= 0:
                tr1_vals = np.hstack((start_zeros,tr1.value))
                tr2_vals = np.hstack((tr2.value, end_zeros))
                times = np.arange(tr2.times.value[0], tr1.times.value[-1]+1/Fs,
                        1/Fs)
                tr1_new = Trace(tr1_vals, times=times, sample_rate=Fs)
                tr2_new = Trace(tr2_vals, times=times, sample_rate=Fs)
            # trace 1 starts before trace 2 and ends after trace 2
            # only add zeros to trace 2
            elif start_diff < 0 and end_diff >= 0:
                tr2_vals = np.hstack((start_zeros, tr2_vals, end_zeros))
                times = tr1.times.value
                tr1_new = np.copy(tr1)
                tr2_new = Trace(tr2_vals, times=times, sample_rate=Fs)
            # trace 1 starts after trace 2 and ends before trace 2
            # swap roles of previous case
            elif start_diff >= 0 and end_diff < 0:
                tr1_vals = np.hstack((start_zeros, tr1_vals, end_zeros))
                times = tr1.times.value
                tr2_new = np.copy(tr2)
                tr1_new = Trace(tr1_vals, times=times, sample_rate=Fs)
            # trace starts before trace 2 and ends before trace 2
            # roles swapped from first case
            elif start_diff < 0 and end_diff < 0:
                tr2_vals = np.hstack((start_zeros,tr2.value))
                tr1_vals = np.hstack((tr1.value, end_zeros))
                times = np.arange(tr1.times.value[0], tr2.times.value[-1],
                        1/Fs)
                tr1_new = Trace(tr1_vals, times=times, sample_rate=Fs)
                tr2_new = Trace(tr2_vals, times=times, sample_rate=Fs)
            dict3[station][channel] = tr1_new + tr2_new

    return dict3

def gaussian_noise2(asd_amp, sample_rate, duration, segdur=None, name=None):
    if segdur is None:
        segdur=duration
    asd = asd_amp * np.ones(100000)
    freqs = np.arange(0,1000, 0.01)
    fi = open('tmp.txt','w')
    for a,f in zip(asd, freqs):
        fi.write('%4.4e %4.4e\n' % (f,a))
    fi.close()
    cmd = 'lalsim-detector-noise -a ./tmp.txt -t %d -d %d -r %d > tmp2.txt' % (duration,
            segdur, sample_rate)
    os.system(cmd)
    final_dat = Trace.read('tmp2.txt', format='txt')

#    os.system('rm tmp.txt tmp2.txt')
    return Trace(final_dat.value, times=final_dat.times,
            sample_rate=final_dat.sample_rate.value, unit=u.m,name=name)

def gaussian_noise(psd_amp, sample_rate, duration, segdur=None, name=None):
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

    return Trace(ts, sample_rate=sample_rate, unit=u.m, name=name)

def get_homestake_channels(type):
    """
    List of all channels saved by Antelope
     HHZ - Channel 1 (z-axis) (100 Hz)
     HHE - Channel 2 (East) (100 Hz)
     HHN - Channel 3 (North) (100 Hz)
     LHZ - Channel 1 (z-axis) (1 Hz)
     LHE - Channel 2 (East) (1 Hz)
     LHN - Channel 3 (North) (1 Hz)
     LCQ - internal clock quality percentage (0-100)
     LCE - clock phase error (microseconds)
     LCC - GPS clock quality
     LPL - clock phase lock loop status
     LCL - time since GPS lock was lost
     QBD - total number of Q330 reboots in last 24 hours
     QBP - logical port buffer percent full from real-time status
     QDG - data gaps (seconds)
     QDL - current data latency (seconds)
     QDR - current total input+output data rate (bits/second)
     QEF - overall communications efficiency (percent)
     QG1 - total number of data gaps in last 1 hour
     QGD - total number of data gaps in last 24 hours
     QID - total number of datalogger ip-address cahnges in last 24 hours
     QLD - total number of comm link cycles in last 24 hours
     QPD - total number of POCs received in last 24 hours
     QRD - total number of bytes read in last 24 hours
     QRT - current run time (seconds)
     QTH - current throttle setting (bits/second0
     QTP - ratio of seconds read to real-time clock
     QWD - total number of bytes written in last 24 hours
     VCO - voltage controlled oscillator value
     VEA - antenna current
     VEC - main system current
     VEP - main system voltage
     VKI - main system temperature
     VM1 - mass position for channel 1
     VM2 - mass position for channel 2
     VM3 - mass position for channel 3
     VPB - percentage of packet buffer full
     VTW - main system opto inputs
     note: all V** channels are updated every 10 seconds
    """
    fast_chans = ['HHE','HHN','HHZ']
    data_chan = ['HHE','HHN','HHZ','LHE','LHN','LHZ']
    status_chan = ['LCE','LCQ','LPL','LCL','LCC','QBD','QBP','QDG','QDL',
            'QDR','QEF','QG1','QGD','QID','QLD','QPD','QRD','QRT',
            'QTH','QTP','QWD','VCO','VEA','VEC','VEP','VKI','VM1',
            'VM2','VM3','VPB','VTW']
    useful_chan = ['HHE','HHN','HHZ','LCE','LCQ','VEP','VKI','VM1','VM2','VM3']

    if (type.lower() == 'all'):
        chan_list = data_chan + status_chan
    elif (type.lower() == 'data'):
        chan_list = data_chan
    elif (type.lower() == 'status'):
        chan_list = status_chan
    elif (type.lower() == 'useful'):
        chan_list = useful_chan
    elif (type.lower() == 'fast_chans'):
        chan_list = fast_chans
    else:
        raise ValueError('type should be \'all\' or \'data\' or \'status\' or \'useful\'.')

    return chan_list

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

