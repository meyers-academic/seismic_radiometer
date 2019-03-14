from gwpy.timeseries import TimeSeries
from numpy.fft import ifft
from gwpy.frequencyseries import FrequencySeries
import astropy.units as u
import numpy as np
from datetime import datetime


def noise_from_psd(length, sample_rate, psd, seed=0, name=None, unit=u.m):
    """ Create noise with a given psd.

    Return noise with a given psd. Note that if unique noise is desired
    a unique seed should be provided.

    Currenlty only works with a single PSD value and assumes white noise.

    Parameters
    ----------
    length : int
        The length of noise to generate in seconds.
    sample_rate : float
       the sample rate of the data
    stride : int
        Length of noise segments in seconds
    psd : FrequencySeries
        The noise weighting to color the noise.
    seed : {0, int}
        The seed to generate the noise.

    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise colored by the given psd.
    """

    if name is None:
        name = 'noise'
    length *= sample_rate
    length = int(length)

    noise_power = psd / 2.
    time = np.arange(length) / sample_rate
    fake_ts = np.random.normal(scale=np.sqrt(noise_power), size=time.shape)

    # THIS WOULD WORK FOR ARBITRARY PSDS ###
    # # fake real and imaginary parts
    # fake_fft_data = 1j*np.random.randn(noise_ts.size) +\
    #     np.random.randn(noise_ts.size)
    # # multply by psd, divide out by amplitude and fix nromalization
    # # so that ifft gives sensible results
    # fake_fft_data *= psd * np.sqrt(2) * np.sqrt(fake_fft_data.size)\
    #     / np.abs(fake_fft_data)
    # newdat = ifft(fake_fft_data)
    return TimeSeries(fake_ts, sample_rate=sample_rate, name=name, unit=unit)


def noise_from_psd2(length, stride, sample_rate, psd, seed=0, name=None,
                    unit=u.m):
    Nsegs = length / stride
    ts = np.zeros(length*sample_rate)
    stride_samps = stride*sample_rate
    for seg in range(Nsegs):
        idx_low = seg*stride_samps
        idx_high = idx_low + stride_samps
        ts[idx_low:idx_high] = noise_from_psd2(stride, sample_rate, psd,
                                               seed=0, name=None).value
    return TimeSeries(ts, name=name, sample_rate=sample_rate, unit=unit)
