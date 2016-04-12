from gwpy.timeseries import TimeSeries
from gwpy.spectrum import Spectrum
import numpy as np
import scipy
from collections import OrderedDict
import glob


class Station(OrderedDict):
    """
    Station class

    initialize station class by:
    >>> frame = 'M-SEISMIC-1125384593-4096.gwf'
    >>> station_name = 'DEAD'
    >>> sta = Station(station_name, st, et)
    >>> print sta
    """

    def __init__(self, st, et, station):
        super(Station, self).__init__()
        self.station = station
        self.st = st
        self.et = et
        self['Z'] = read_frame('M-SEISMIC-1125384593-4096.gwf',
                               station + ':HHZ', st=st, et=et)
        self['N'] = read_frame('M-SEISMIC-1125384593-4096.gwf',
                               station + ':HHN', st=st, et=et)
        self['E'] = read_frame('M-SEISMIC-1125384593-4096.gwf',
                               station + ':HHE', st=st, et=et)
        self.location = self['Z'].location


class StationArray(OrderedDict):
    """
    StationArray class

    >>> frame = 'M-SEISMIC-1125384593-4096.gwf'
    >>> station_names = ['DEAD','D4850']
    >>> arr = StationArray(station_names, st, et)
    >>> print arr
    >>> arr[station_name[0]]['Z']
    """

    def __init__(self, st, et, stations):
        super(StationArray, self).__init__()
        self.st = st
        self.et = et
        self.stations = stations
        for station in stations:
            self[station] = Station(st, et, station)

    def coherence(self, fftlength=None, window='hanning',
                  **kwargs):
        COH = OrderedDict()
        for ii in range(len(self.keys())):
            for jj in range(ii, len(self.keys())):
                if ii == jj:
                    continue
                newkey = self.keys()[ii] + '-' + self.keys()[jj]
                key1 = self.keys()[ii]
                key2 = self.keys()[jj]
                COH[newkey] = self[key1]['Z'].coherence(self[key2]['Z'],
                                                        fftlength=fftlength,
                                                        window=window,
                                                        stacktype='ts',
                                                        **kwargs)
        return COH


class Trace(TimeSeries):
    """class for doing seismic data analysis, inherited from gwpy TimeSeries"""

    def hilbert(self):
        """
        Performs hilbert transform to get envelope of TS data.

        Parameters
        ----------

        Returns
        -------
        TS : `Trace`
            envelope of data. modulus of output of
            scipy.signal.hilbert, which is actual full analytic
            extension of data (not just what woudl normally
            be considered hilbert transform)
        """
        # hilbert transform
        arr = scipy.signal.hilbert(self.value)

        # class assignment stuff...
        TS = Trace(arr)
        TS.__dict__ = self.copy_metadata()
        return TS.detrend()

    def smooth(self, width=1):
        """
        Smooths data by convolving data with ones...
        still not sure how exactly this works and seems to
        produce some edge effects

        Parameters
        ----------
        width : `int`, optional, default=1,
            number of seconds or Hz to use in convolution.

        Returns
        -------
        smoothed : `Trace`
            Smoothed time series trace.
        """

        TS = np.abs(self)

        # turn width into # of samples
        width = width * TS.sample_rate.value

        # get window
        window = np.ones((2 * width))

        # do convolution
        TS = np.convolve(TS, window / (2 * width), 'same')
        TS = Trace(TS)
        TS.__dict__ = self.copy_metadata()
        return TS.detrend()

    def renormalization(self, Ns=None, type='water_level'):
        """
        Does renormalization to get rid of things like
        EQs.

        'weighted_renorm': Calculates weights based on
        earthquake-band bandpass filter,
        applies those weights to raw data. (requires
        number of seconds with which to calculate weights)

        'water_level': calculates envelope of data using smoothed
        modulus of hilbert transform and normalizes by that
        smoothed envelope

        'bit': one bit normalization. return sign of detrended
        data

        Parameters:
        -----------
        Ns : int, optional
            Number of seconds to include in
            renormalization. Optimal is half of max period
            of planned bandpass. Must be set if type is
            'weighted_renorm'.
        Type : str, optional, default = 'water_level'
            Type of renormalization. Default
            is water level renormalization.
            Other options: 'bit', 'weighted_renorm'

        Returns:
        --------
            normed : `Trace` object, weighted and renormalized
        """

        # weighted renormalization
        if type == 'weighted_renorm':
    # need width to use for weights
            if Ns is None:
                raise ValueError(
                    'if type is weighted_renorm you\
                    need number of seconds to renormalize by!')
            TS = self
            # apply EQ bandpass filter
            EQ_TS = TS.bandpass(0.03, 0.1)
            Nsamps = Ns * EQ_TS.sample_rate.value
            idx = 0
            normed = np.zeros(EQ_TS.size)
            # get weights from bandpassed data, apply them to raw data
            for datum in TS.value:
                if idx >= Nsamps and idx + Nsamps <= EQ_TS.value.size:
                    normed[idx] = datum / \
                        np.mean(np.abs(EQ_TS.value[idx - Nsamps:idx + Nsamps]))
                elif idx - Nsamps < 0:
                    normed[idx] = datum / \
                        np.mean(np.abs(EQ_TS.value[:idx + Nsamps]))
                elif idx + Nsamps > EQ_TS.value.size:
                    normed[idx] = datum / \
                        np.mean(np.abs(EQ_TS.value[idx - Nsamps:]))
                idx += 1
            normed = Trace(
                normed, sample_rate=TS.sample_rate,
                name=TS.name)
            return normed.detrend()
        # water level renormalization
        elif type == 'water_level':
            # take hilbert transform
            hil = self.hilbert()

            # get envelope
            env = np.abs(hil)

            # smooth envelope (take one sample each second)
            env = env.smooth()

            # normalize by envelope
            TS = self / env
            return TS.detrend()

        elif type == 'bit':
            TS = np.sign(self.detrend())
            TS = Trace(TS)
            TS.__dict__ = self.copy_metadata()
            return TS.detrend()

    def fft_new(self, **kwargs):
        """
        Calculates fft (keeps negative frequencies as well...
        we need them for ambient noise cross correlation).

        NOTE: This is renormalized to be the correct spectrum,
        however that means that you cannot just use
        numpy.fft.ifft(self.fft_new(window=None))
        to get back the original timeseries.

        >>> data1 = read_frame(frame, channel)
        >>> TS_old = np.fft.ifft(data1.size * data1.fft_new(window=None))
        >>> data1 == TS_old

        self.fft() uses the same normalization, but does not offer
        whitening and windowing like this function does.

        can do whitening if you want

        Parameters:
        -----------
        whiten: `bool`, optional
            Whitens spectrum on both sides.

        Returns:
        --------
        fft : `Spec`, fft
            Whitened if wanted
        """
        try:
            kwargs['whiten']
        except KeyError:
            kwargs['whiten'] = False
        try:
            kwargs['window']
        except KeyError:
            kwargs['window'] = 'hanning'

        # whiten
        # if kwargs['whiten']:
        #     Nsecs = self.value.size / self.sample_rate.value
        #     TS = self.whiten(1. / 8 * Nsecs, 1. / 16 * Nsecs)
        # else:
        TS = self
        # window and fft...
        if kwargs['window'] == 'hanning':
            window = np.hanning(TS.value.size)
            fft = np.fft.fft(
                TS.value * window) / (TS.size)
            freqs = np.fft.fftfreq(
                TS.value.size, d=(1. / TS.sample_rate.value))
        else:
            fft = np.fft.fft(
                TS.value) / TS.size
            freqs = np.fft.fftfreq(
                TS.value.size, d=(1. / TS.sample_rate.value))
        # put things in the right order
        fft = Spec(
            fft, f0=freqs[0], df=(freqs[1] - freqs[0]),
            name=self.name, epoch=self.epoch)
        if kwargs['whiten']:
            fft = fft.whiten(width=1)

        return fft

    def coherence_calc(self, tr, whiten=False, bandpass=None, flow=1e-4,
                       fhigh=50, normtype=None, normlen=None, window='hanning',
                       fftlength=None, overlap=None, outtype='ts'):
        """
        Calculates coherence between two traces. Will do spectral whitening,
        normalize it.

        Parameters:
        -----------
        tr : `Trace`
            Trace time series to calculate coherence with.

        Return:
        -------
        ifft : `Trace`
            Trace time series of coherence
        """
        if bandpass is None and normtype is None:
            new1 = self.fft_new(whiten=whiten, window=window)
            new2 = tr.fft_new(whiten=whiten, window=window)
        elif bandpass is None and normtype is not None:
            new1 = self.renormalization(
                Ns=normlen, type=normtype).fft_new(whiten=whiten,
                                                   window=window)
            new2 = tr.renormalization(
                Ns=normlen, type=normtype).fft_new(whiten=whiten,
                                                   window=window)
        elif normtype is not None and bandpass is not None:
            new1 = self.bandpass(flow, fhigh).renormalization(
                Ns=normlen, type=normtype).fft_new(whiten=whiten,
                                                   window=window)
            new2 = tr.bandpass(flow, fhigh).renormalization(
                Ns=normlen, type=normtype).fft_new(whiten=whiten,
                                                   window=window)
        elif normtype is None and bandpass is True:
            new1 = self.bandpass(flow, fhigh).fft_new(
                whiten=whiten, window=window)
            new2 = tr.bandpass(flow, fhigh).fft_new(
                whiten=whiten, window=window)

        coh = np.conj(new1) * new2 / (np.abs(new1) * np.abs(new2))
        # do ifft, but need to rearrange frequencies
        coh_ts = np.fft.ifft(coh)
        coh_ts = TimeSeries(coh_ts, name='coherence TS between %s and %s' % (
            self.channel, tr.channel), sample_rate=(self.sample_rate))
        if outtype == 'ts':
            return coh_ts
        if outtype == 'spec':
            return np.fft.fftshift(coh)
        if outtype == 'components':
            return new1, new2

    def coherence(self, tr, window='hanning', fftlength=None, stacktype='freq',
                  **kwargs):
        """
        Calculate coherence between self and `tr` trace object.

        Parameters
        ----------
        tr : `Trace`
            object with with to calculate coherence
        fftlength : `int`
            length of ffts to take when calculating coherence
        stacktype : `str`
            method used to stack coherences. options are 'ts' and 'freq'
            'ts' averages resultant ifft'ed timeseries together. 'freq'
            averages csds and psds individually and takes ratio of them at
            the end and then takes ifft.

        Returns
        -------
        coh_ts : `Trace`
            coherence timeseries with acausal followed by causal times
        """
        if fftlength is None:
            fftlength = self.size
        Nsamps = fftlength * self.sample_rate.value

        # stack
        if window == 'hanning' and fftlength <= 2 * self.size:
            nsteps = 2 * self.size / Nsamps
        else:
            nsteps = self.size / Nsamps
        for step in range(int(nsteps) - 1):
            idx1 = step * Nsamps / 2.
            idx2 = idx1 + Nsamps
            if stacktype is 'freq':
                new1, new2 = self[idx1:idx2].coherence_calc(
                    tr[idx1:idx2], outtype='components', window=window,
                    **kwargs)
                if step == 0:
                    csd = np.conj(new1) * new2
                    asd1 = np.abs(new1)
                    asd2 = np.abs(new2)
                else:
                    csd = (csd * step + new1 * new2) / (step + 1)
                    asd1 = np.sqrt(
                        (asd1 ** 2 * step + np.abs(new1)**2) / (step + 1))
                    asd2 = np.sqrt(
                        (asd2 ** 2 * step + np.abs(new2)**2) / (step + 1))
            elif stacktype is 'ts':
                if step == 0:
                    coh_ts = self[idx1:idx2].coherence_calc(
                        tr[idx1:idx2], outtype='ts', window=window,
                        **kwargs)
                else:
                    coh_ts_temp = self[idx1:idx2].coherence_calc(
                        tr[idx1:idx2], outtype='ts', window=window,
                        **kwargs)
                    coh_ts = (coh_ts * step + coh_ts_temp) / (step + 1)

        if stacktype is 'freq':
            coh = csd / (asd1 * asd2)
            coh_ts = np.fft.ifft(coh)
        deltaXvec = self.location - tr.location
        N = coh_ts.size
        coh_ts = np.hstack((coh_ts[N / 2:], coh_ts[:N / 2]))
        coh_ts = Trace(coh_ts, name='coherence TS between %s and %s' % (
            self.channel, tr.channel), sample_rate=(self.sample_rate))
        coh_ts.deltax = deltaXvec
        coh_ts.x0 = -coh_ts.times[-1] / 2.
        return coh_ts.real

    def get_location(self):
        """
        Gets location of a specific station based on channel name associated
        with trace.

        Parameters
        ----------
        none

        Returns
        -------
        location : `numpy array`
            [x,y,z] location of station based on channel name
        """

        xyz_list = {'DEAD': [599316.496208, 4915135.795515, 1498],
                    'LHS': [597654.698101, 4911177.756239, 1684],
                    'ORO': [599454.495625, 4910782.727115, 1543],
                    'ROSS': [599029.335575, 4910954.029719, 1628],
                    'RRDG': [598383.512647, 4912544.118885, 1677],
                    'SHL': [602888.695761, 4907880.584008, 1772],
                    'TPK': [595816.121345, 4910428.385396, 1740],
                    'WTP': [600233.873095, 4911950.092981, 1555],
                    'YATES': [599503.542430, 4911750.052702, 1625],
                    '300': [599082.941268, 4911099.273511, 1505.1],
                    '800': [598921.845912, 4911207.932696, 1350],
                    '1700': [598954.462337, 4911686.159936, 1073.8],
                    'A2000': [598644.641050, 4911614.812799, 983.3],
                    'B2000': [598791.293544, 4911405.938094, 983.3],
                    'C2000': [598532.222831, 4911668.666166, 983.1],
                    'D2000': [598114.948236, 4911851.254988, 983],
                    'E2000': [597894.603780, 4912192.359339, 983.1],
                    'A4100': [599356.477888, 4910936.776895, 342.5],
                    'C4100': [599568.994798, 4911639.949104, 342.3],
                    'D4100': [599549.979057, 4910795.291477, 342.4],
                    'A4850': [598646.478984, 4910437.175145, 115.2],
                    'B4850': [598987.461619, 4911086.715912, 114.9],
                    'C4850': [599409.907522, 4911093.130999, 114.6],
                    'D4850': [599581.886292, 4911840.127688, 115.2]}
        staname = self.channel.name.split(':')[0]
        return np.asarray(xyz_list[staname])


class Spec(Spectrum):
    def whiten(self, width=1):
        """
        whitens spectrum by getting a smoothed version of the
        absolute value of the spectrum and dividing out by that smoothed
        version.

        Parameters
        ----------
        width : `int`, kwarg, optional, default=1
            width over which to smooth (in Hz)

        Returns
        -------
        out : `Spec`, whitened spectrum
        """

        S = self
        env = S.smooth(width=width)
        return S / env

    def smooth(self, width=1):
        """
        Smooths spectrum by convolving data with ones of
        proper length. Averages over width

        Parameters
        ----------
        width : `int`, optional, default=1,
            Number of seconds or Hz to use in convolution.

        Returns
        -------
        smoothed : `Trace`
            Smoothed time series trace.
        """

        S = np.abs(self)

        # turn width into # of samples
        width = width * (1 / S.df.value)

        # get window
        window = np.ones((2 * width))

        # do convolution
        S = np.convolve(S, window / (2 * width), 'same')
        S = Spec(S)
        S.__dict__ = self.copy_metadata()
        return S


def fetch(st, et, channel, framedir='./'):
    """
    fetch data based on location of frames

    Parameters
    ----------
    st : `int`
        start time (GPS time)
    et : `int`
        end time (GPS time)
    channel : `channel`
        channel to load data for

    Returns
    -------
    TS : `Trace`
        Trace object containing data between start and end times
    """
    # uncomment when not testing
    # for looping over directories where frames
    # are located
    # st_dir = int(str(st)[:5])
    # et_dir = int(str(et)[:5])
    # dirs = np.arange(st_dir, et_dir + 1)
    # for dir in dirs:
    print 'FRAME READING TESTING MODE!!!'
    files = sorted(glob.glob(framedir + '*.gwf'))
    vals = np.asarray([])

    for file in files:
        fst = file.split('-')[2]
        dur = file.split('-')[-1][:-4]
        # start is before frame start, end is before frame end
        # we want to load fst -> et
        if st <= fst and et <= fst + dur:
            val = read_frame(file, channel, st=fst, et=et)
            vals = np.hstack((vals, val.value))
        # start is after frame start, end is before frame end
        # we want to load only st -> et
        elif st >= fst and et <= fst + dur:
            val = read_frame(file, channel, st=st, et=et)
            vals = np.hstack((vals, val.value))
        # start is after frame start end is after or equal to frame end
        # we want to load st -> fst + dur
        elif st >= fst and et >= fst + dur:
            val = read_frame(file, channel, st=st, et=fst + dur)
        # start is before frame start, end is after frame end
        # load fst -> fst + dur (whole frame)
        elif st <= fst and et >= fst + dur:
            val = read_frame(file, channel, st=fst, et=fst + dur)
            vals = np.hstack(vals, val.value)
        TS = Trace(vals, x0=st, dx=val.dx, name=val.name, channel=val.channel)
        loc = TS.get_location()
        TS.location = loc
        return TS


def read_frame(frame, channel, st=None, et=None, cfac=1.589459e-9):
    """
    reads ligo frames

    Parameters
    ----------
    frame : `str`
        filepath to a frame
    channel : `str`
        channel in the frame to load
    st : `int`, date string, optional
        optional start time. defaults to beginning
        of frame
    et : `int ,date string, optional
        optional end time. defaults to end
        of frame

    Returns
    -------
    TS : `Trace`
        time series trace
    """

    if st is not None and et is not None:
        d1 = cfac * Trace.read(frame, channel, st, et).detrend()
        d2 = TimeSeries.read(frame, channel, st, et)
    else:
        d1 = cfac * Trace.read(frame, channel).detrend()
        d2 = TimeSeries.read(frame, channel)
    d1.__dict__ = d2.copy_metadata()
    d1.location = d1.get_location()
    return d1
