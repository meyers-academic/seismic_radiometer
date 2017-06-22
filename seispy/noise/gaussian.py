import lal
from lalsimulation import SimNoise
from gwpy.timeseries import TimeSeries
from gwpy.frequencyseries import FrequencySeries
import astropy.units as u
import numpy as np
from datetime import datetime


def noise_from_psd(length, sample_rate, psd, seed=0, name=None, unit=u.m):
    """ Create noise with a given psd.

    Return noise with a given psd. Note that if unique noise is desired
    a unique seed should be provided.
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
        name='noise'
    length *=sample_rate
    length = int(length)

    noise_ts = TimeSeries(np.zeros(int(length)),
            sample_rate=sample_rate, name=name, unit=unit)
    if seed == 0 or seed is None:
        now = datetime.now()
        seed = now.year+now.day+now.hour+now.minute+now.second+now.microsecond

    randomness = lal.gsl_rng("ranlux", seed)

    N = int (sample_rate / psd.df.value)
    n = N/2+1
    stride = N/2

    if n > len(psd):
        raise ValueError("PSD not compatible with requested delta_t")
    psd = (psd[0:n]).to_lal()
    psd.data.data[n-1] = 0
    segment = TimeSeries(np.zeros(N), sample_rate=sample_rate).to_lal()
    length_generated = 0
    newdat=[]

    SimNoise(segment, 0, psd, randomness)
    while (length_generated < length):
        if (length_generated + stride) < length:
            newdat.extend(segment.data.data[:stride])
        else:
            newdat.extend(segment.data.data[0:length-length_generated])

        length_generated += stride
        SimNoise(segment,stride, psd, randomness)
    return TimeSeries(newdat, sample_rate=sample_rate, name=name, unit=unit)

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

