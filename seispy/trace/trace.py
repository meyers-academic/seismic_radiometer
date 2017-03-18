from gwpy.timeseries import TimeSeries
import numpy as np
import scipy
import glob
from gwpy.time import *
from astropy import units as u


class Trace(TimeSeries):
    """class for doing seismic data analysis, inherited from gwpy TimeSeries"""

    def to_obspy(self):
        import obspy.core.trace as tr
        """
        Take this trace to an obspy trace

        Returns
        -------
        """
        stats = tr.Stats()
        stats['sampling_rate'] = self.sample_rate.value
        stats['location'] = self.get_location()
        stats['starttime'] = self.times.value[0]
        stats['npts'] = self.times.value.size
        stats['channel'] = self.channel.name.split(':')[1]
        stats['station'] = self.channel.name.split(':')[0]
        return tr.Trace(self.value, stats)

    def hilbert(self):
        """
        TODO: need a unittest

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
        TS = Trace(arr, name=self.name, channel=self.channel,
                   sample_rate=self.sample_rate, unit=self.unit,
                   epoch=self.epoch)
        return TS.detrend()

    def taper(self, **kwargs):
        """
        TODO: unittest

        Returns a new (deep) copy of the object with a sin^2 taper
        applied to both ends of the data.

        Takes kwargs of 'v_min, v_max, dist' or a
        :class:`seispy.event.Event` class
        """

        # Make copy.
        arr = np.copy(self.value)
        dt = 1 / self.sample_rate.value

        # Set up taper lengths
        if all(k in kwargs for k in ('v_min','v_max','dist')):
            dur_t1 = float(kwargs['dist']/kwargs['v_max'])
            dur_t2 = (arr.size*dt) - (float(kwargs['dist']/kwargs['v_min']) + 60)
            # Handle case where arrival time is after the end of the data.
            # Taper just last 10 seconds.
            if (dur_t2 <= 10):
                dur_t2 = 10
        elif ('event' in kwargs):
            dur_t1 = kwargs['event'].taper_start
            dur_t2 = kwargs['event'].taper_end
        else:
            raise ValueError('kwargs should be v_min, v_max, dist OR event=seisEv.')


        # Make beginning taper.
        nSamp_t1 = int(np.ceil(dur_t1/dt))
        samples1 = np.arange(0,nSamp_t1,step=1)
        taper1 = np.sin(2*np.pi*samples1/(4*(nSamp_t1-1)))

        # Make end taper.
        # Check if end arrival time is aft
        nSamp_t2 = int(np.ceil(float(dur_t2)/dt))
        samples2 = np.arange((nSamp_t2-1), -1,step=-1)
        taper2 = np.sin(2*np.pi*samples2/(4*(nSamp_t2-1)))

        # Taper data.
        times = self.times.value
        arr[0:nSamp_t1] *= taper1
        arr[-nSamp_t2:] *= taper2
        newtimes = np.roll(times, -nSamp_t2)[:(nSamp_t1+nSamp_t2)]
        newtrace = Trace(arr, times=newtimes, sample_rate=self.sample_rate,
                         channel=self.channel, name=self.name)
        newtrace.tapered = True

        return newtrace

    def smooth(self, width=1):
        """
        TODO: need a unittest

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
        arr = np.convolve(TS.value, window / (2 * width), 'same')
        TS = Trace(arr, name=self.name, channel=self.channel,
                   sample_rate=self.sample_rate, unit=self.unit,
                   epoch=self.epoch)
        return TS.detrend()

    def gaussian_filter(self, filt_freq, filt_width):
        """
        Parameters
        ----------
        filt_freq : `float`
            frequency of center of our filter
        filt_width : `float`
            width of our gaussian filter in Hz

        Returns
        -------
        new_ts : `seispy.trace.Trace`
            new trace with filtered data
        """
        new_fft = self.fft()
        # create a gaussian filter
        exp_array = -(filt_freq - new_fft.frequencies.value) ** 2 / ((filt_width * filt_freq) ** 2)
        # multiply fft by our gaussian filter
        new_fft *= np.exp(exp_array)
        # take ifft
        new_ts = new_fft.ifft()
        # set up and return new Trace object
        new_trace = Trace(new_ts.value, sample_rate=self.sample_rate,
                          name=self.name, epoch=self.epoch,
                          unit=self.unit)
        return new_trace

    def vel2disp(self):
        """
        Change this trace from velocity to displacement.

        Returns
        -------
        new_trace : `seispy.trace.Trace`
            new trace with units of input units * seconds
        """
        # get fft
        new_fft = self.fft()
        # divide by frequencies
        freqs = new_fft.frequencies.value
        freqs[0] = 1
        new_fft *= (1 / new_fft.frequencies.value)
        # take ifft
        new_ts = new_fft.ifft()
        # get new Trace object
        new_trace = Trace(new_ts.value, sample_rate=self.sample_rate, name=self.name, epoch=self.epoch, unit=self.unit*u.s)
        return new_trace


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
#            TS.__dict__ = self.copy_metadata()
            return TS.detrend()

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
    st_dir = int(str(st)[:5])
    et_dir = int(str(et)[:5])
    dirs = np.arange(st_dir, et_dir + 1)
    files = []
    for directory in dirs:
        # print 'FRAME READING TESTING MODE!!!'
        loaddir = '%s/M-%d/' % (framedir, directory)
        new_files = sorted(glob.glob(loaddir + '/*.gwf'))
        files.extend(new_files)
    vals = np.asarray([])
    if len(files)==0:
        raise ValueError('No files found...we looked here: %s' % loaddir)

#    for file in files:
    for ii in range(len(files)):
        file = files[ii]
        # start is before frame start, end is before frame end
        # we want to load fst -> et
        fst = int(file.split('-')[-2])
        dur = int(file.split('-')[-1][:-4])
        if st <= fst and et <= fst + dur and et >= fst:
            val = read_frame(file, channel, st=fst, et=et)
            vals = np.hstack((vals, val.value))
        # start is after frame start, end is before frame end
        # we want to load only st -> et
        elif st >= fst and et <= fst + dur:
            val = read_frame(file, channel, st=st, et=et)
            vals = np.hstack((vals, val.value))
        # start is after frame start end is after or equal to frame end
        # we want to load st -> fst + dur
        elif st >= fst and st < (fst + dur) and et >= fst + dur:
            val = read_frame(file, channel, st=st, et=fst + dur)
            vals = np.hstack((vals, val.value))
        # start is before frame start, end is after frame end
        # load fst -> fst + dur (whole frame)
        elif st <= fst and et >= fst + dur:
            val = read_frame(file, channel)
            vals = np.hstack((vals, val.value))
        else:
            continue
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
    else:
        d1 = cfac * Trace.read(frame, channel).detrend()
    d1.location = d1.get_location()
    return d1


def read_mseed(f, starttime=None, endtime=None):
    """
    reads miniseed file between start time and end time.

    NOTE: No error is given if times are outside of range
    of data in file. Data in file is clipped to start and
    end times and if data in file is subset of specified
    interval then all data is returned. This is because
    this function is used primarily as a tool for `fetch_mseed()`
    and this makes life much easier.

    Parameters
    ----------
    f : `str`
        miniseed file to load
    starttime : `int`, optional
        GPS time. Start time. defaults
        to start time of file to load
    endtime : `int`, optional
        GPS time. End time. Defaults to
        end time of file to load.
    """
    from obspy import (read, UTCDateTime)
    from gwpy.time import tconvert

    # if no start time or end time, read the whole
    # damn thing
    data = read(f)
    data_st = float(tconvert(UTCDateTime(data[0].stats.starttime))) - 315964783
    data_et = float(tconvert(UTCDateTime(data[0].stats.endtime))) - 315964783
    dx = 1. / data[0].stats.sampling_rate
    if starttime is None and endtime is None:
        print('No start times specified, reading whole file...')
    elif isinstance(starttime, int) and isinstance(endtime, int):
        # check start/end times match file we're reading...useful
        # for fetching
        if starttime <= data_st:
            st = UTCDateTime(tconvert(data_st))
        else:
            st = UTCDateTime(tconvert(starttime))
        if endtime - dx >= data_et:
            et = UTCDateTime(tconvert(data_et - dx))
        else:
            et = UTCDateTime(tconvert(endtime - dx))
        data = read(f, starttime=st, endtime=et)

    tr = data[0]
    trace = Trace(tr.data, x0=starttime, dx=1. / tr.stats.sampling_rate,
                  channel='%s:%s' % (data[0].stats.station, data[0].stats.channel))
    return trace


def find_mseed_files(channel, st, et, basedir='./'):
    if isinstance(st, int):
        st = from_gps(st)
    if isinstance(et, int):
        et = from_gps(et)
    day1 = st.timetuple().tm_yday
    dayend = et.timetuple().tm_yday
    files = []
    station = channel.split(':')[0]
    chan = channel.split(':')[1]
    # check if we change years...
    if st.year < et.year:
        if et.year - st.year == 1:
            for ii in range(day1, 366):
                files.append('%d/%s/%03d/%s.X6..%s.%d.%03d' %
                             (st.year, basedir, ii, station, chan, st.year, ii))
            for jj in range(1, dayend + 1):
                files.append('%d/%s/%03d/%s.X6..%s.%d.%03d' %
                             (et.year, basedir, jj, station, chan, et.year, jj))
        else:
            raise ValueError('for now cant go across multiple years...')
    else:
        for ii in range(day1, dayend + 1):
            files.append('%d/%s/%03d/%s.X6..%s.%d.%03d' %
                         (st.year, basedir, ii, station, chan, st.year, ii))
    return files


def fetch_mseed(channel, st, et, basedir='./', cfac=1.589459e-9):
    """
    fetch miniseed data from database with toplevel directory
    of `basedir`.

    Paramters
    ---------
    channel : `str`
        Station/channel pair. e.g. 'D4850:HHZ' to match frame
        fetching method
    st : `int`
        LIGO gps time for start
    et : `int`
        LIGO gps time for end (will not include sample starting on end time)
    basedir `str`, optional, default='./'
        Base directory for database of antelope data
    cfac : `float`
        Calibration constant factor to multiply by data from minised files.

    Returns
    -------
    data : :Trace:
        Detrended, calibrated trace object
    """
    files = find_mseed_files(channel, st, et, basedir=basedir)
    data = np.array([])
    st = st
    # test open to get sample rate of channel
    dat_test = read_mseed(files[0], st, et)
    if from_gps(et).timetuple().tm_yday - from_gps(et - dat_test.dx.value).timetuple().tm_yday > 0:
        files = files[:-1]
    for f in files:
        dat = read_mseed(f, st, et)
        data = np.hstack((data, dat.value))
    data = cfac * Trace(data, x0=st, dx=dat.dx, channel=channel)
    return data.detrend()
