from __future__ import division, print_function
from collections import OrderedDict
from ..utils import *
import astropy.units as u
from ..noise import gaussian
from ..trace import Trace, fetch
import numpy as np
from gwpy.frequencyseries import FrequencySeries
from .station import StationArray, homestake
import healpy as hp


class Seismometer(OrderedDict):
    """
    Station data
    """

    def __copy__(self):
        new = type(self)()
        for key in self.keys():
            new[key] = self[key].copy()

    @classmethod
    def fetch_data(cls, station_name, st, et, framedir='./',
                   chans_type='useful'):
        """
        fetch data for this seismometer

        Parameters
        ----------
        station_name : `str`
            name of station for which we're fetching data
        st : `int`
            start time
        et : `int`
            end time
        framedir : TODO, optional
            frame base directory
        chans_type : `str`
            type of channels to return. 'fast_chans`
            is for our HHZ, HHE, HHN chans and is a
            safe bet. `useful_chans` are for if we're
            interested in looking at clock quality

        Returns
        -------
        seismometer : :class:`stationdata.Seismometer`
            seismometer data
        """
        seismometer = cls()
        chans = get_homestake_channels(chans_type)
        for chan in chans:
            seismometer[chan] = fetch(st, et, station_name + ':' + chan,
                                      framedir=framedir)
        return seismometer

    @classmethod
    def initialize_all_good(cls, duration, chans_type='useful', start_time=0,
                            name='seismometer', location=None, Fs=None):
        """

        Parameters
        ----------
        duration : `int`
            length of data we're using
        chans_type : `str`
            type of channels
        start_time : `int`
            start time of initialized data
        name : `str`
            name for this data
        location : `numpy.ndarray`
            location array

        Returns
        -------
        seismometer : :class:`seispy.station.stationdata.Seismometer`
            seismometer data with zeros initialized everywhere

        Initialize a seismometer with all good
        information in all channels and zeros in
        the data channels.
        channel | srate
        HHN 100
        HHE 100
        LCE 1
        HHZ 100
        LCQ 1
        VKI 0.1000000
        VEP 0.1000000
        VM1 0.1000000
        VM3 0.1000000
        VM2 0.1000000
        """
        if location is None:
            location = [0, 0, 0]
        name = str(name)
        if Fs is None:
            Fs = 100
        seismometer = Seismometer()
        seismometer['HHE'] = Trace(np.zeros(int(duration * Fs)),
                                   sample_rate=Fs * u.Hz,
                                   epoch=start_time, name=name + ' East',
                                   unit=u.m, channel='%s:%s' % (name, 'HHE'))
        seismometer['HHN'] = Trace(np.zeros(int(duration * Fs)),
                                   sample_rate=Fs * u.Hz,
                                   epoch=start_time, name=name + ' North',
                                   unit=u.m, channel='%s:%s' % (name, 'HHN'))
        seismometer['HHZ'] = Trace(np.zeros(int(duration * Fs)),
                                   sample_rate=Fs * u.Hz,
                                   epoch=start_time, name=name + ' Vertical',
                                   unit=u.m, channel='%s:%s' % (name, 'HHZ'))
        if chans_type == 'fast_chans':
            for chan in seismometer.keys():
                seismometer[chan].location = location
            return seismometer
        else:
            seismometer['LCQ'] = Trace(100 * np.ones(int(duration * 1)),
                                       sample_rate=1 * u.Hz, epoch=start_time,
                                       name=name + ' Clock Quality',
                                       channel='%s:%s' % (name, 'LCQ'))
            seismometer['LCE'] = Trace(np.zeros(int(duration * 1)),
                                       sample_rate=1 * u.Hz, epoch=start_time,
                                       name=name + ' Clock Phase\
                    Error', channel='%s:%s' % (name, 'LCE'))
            seismometer['VM1'] = Trace(np.zeros(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz,
                                       epoch=start_time, name=name + ' Mass\
                    Position Channel 1', channel='%s:%s' % (name, 'VM1'))
            seismometer['VM2'] = Trace(np.zeros(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz,
                                       epoch=start_time, name=name + ' Mass\
                    Position Channel 2', channel='%s:%s' % (name, 'VM2'))
            seismometer['VM3'] = Trace(np.zeros(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz,
                                       epoch=start_time, name=name + ' Mass\
                    Position Channel 3', channel='%s:%s' % (name, 'VM3'))
            seismometer['VEP'] = Trace(13 * np.ones(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz,
                                       epoch=start_time, name=name + ' System\
                    Voltage', channel='%s:%s' % (name, 'VEP'))
            seismometer['VKI'] = Trace(np.zeros(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz,
                                       epoch=start_time, name=name + ' System\
                    temperature', channel='%s:%s' % (name, 'VKI'))
            # set location
            for chan in seismometer.keys():
                seismometer[chan].location = location
        return seismometer

    def to_obspy(self):
        import obspy
        mystream = obspy.core.stream.Stream()
        for key in self.keys():
            mystream += self[key].to_obspy()
        return mystream

    def rotate_RT(self, bearing):
        """

        Parameters
        ----------
        bearing : `float`
            degrees from north

        Returns
        -------
        seismometer : `seispy.station.stationdata.Seismometer`
            seismometer with new channels for transverse and radial related
            to a single source
        """
        theta = np.radians(bearing)
        dataT = self['HHN'].value * np.cos(theta) -\
            self['HHN'].value * np.sin(theta)
        dataR = self['HHE'].value * np.sin(theta) +\
            self['HHE'].value * np.cos(theta)
        sample_rate = self['HHZ'].sample_rate
        epoch = self['HHZ'].epoch
        self['HHT'] = Trace(dataT, sample_rate=sample_rate, epoch=epoch,
                            name='Transverse relative to %4.2f' % bearing,
                            channel=self['HHN'].channel)
        self['HHR'] = Trace(dataR, sample_rate=sample_rate, epoch=epoch,
                            name='Radial relative to %4.2f' % bearing,
                            channel=self['HHN'].channel)


class SeismometerArray(OrderedDict):
    """
    Seismometer array object

    """

    # TODO add taper, filter, etc.
    def to_obspy(self):
        from obspy.core.stream import Stream
        mystream = Stream()
        for key in self.keys():
            mystream += self[key].to_obspy()

        return mystream

    @classmethod
    def fetch_data(cls, st, et, framedir='./', chans_type='useful',
                   stations=None):
        """
        Parameters
        ----------
        st : `int`
            start time
        et : `int`
            end time
        framedir : `string`, optional
            top level frame directory
        chans_type : `type of chans to load`, optional

        Returns
        -------
        seismometer_array : :class:`seispy.station.stationdata.SeismometerArray`
            seismometer array
        """
        # TODO: Docstring for fetch_data.
        if stations is not None:
            arr = cls.initialize_all_good(homestake(stations), et - st,
                                          chans_type=chans_type, start_time=st)
        else:
            arr = cls.initialize_all_good(homestake(), et - st,
                                          chans_type=chans_type, start_time=st)

        for station in arr.keys():
            arr[station] = Seismometer.fetch_data(station, st, et,
                                                  framedir=framedir,
                                                  chans_type=chans_type)
        return arr

    def rotate_RT(self, bearing):
        for key in self.keys():
            self[key].rotate_RT(bearing)

    def filter_by_depth(self, depth):
        new_dict = self.copy()
        for key in new_dict.keys():
            if new_dict[key]['HHE'].get_coordinates()[-1] != depth:
                _ = new_dict.pop(key)
        return new_dict

    def filter_by_channels(self, channels):
        newdict = self.copy()
        for key in newdict.keys():
            for key2 in newdict[key].keys():
                if not any(key2 in c for c in channels):
                    _ = newdict[key2].pop(key2)
        return newdict

    def mean_rz_phase(self, flow, fhigh, *args, **kwargs):
        phases = SeismometerArraySpectrum()
        phases_stds = SeismometerArraySpectrum()
        for sta in self.keys():
            phases[sta] = OrderedDict()
            phases_stds[sta] = OrderedDict()
            try:
                phase = self[sta]['HHR'].csd_spectrogram(self[sta]['HHZ'],
                                                         *args,
                                                         **kwargs)
            except KeyError:
                raise KeyError('No radial channel defined yet')
            freqs = phase.frequencies.value
            phase = np.angle(phase)
            phase_avg = FrequencySeries(np.mean((phase), axis=0),
                                        frequencies=freqs)
            phase_stds = FrequencySeries(np.std((phase), axis=0),
                                        frequencies=freqs)
            phases[sta]['phase'] = phase_avg
            phases_stds[sta]['phase'] = phase_stds
        return phases, phases_stds

    def spectrogram(self, *args, **kwargs):
        spec_dict = SeismometerArraySpecgram()
        for loc in self.keys():
            spec_dict[loc] = {}
            for chan in self[loc].keys():
                spec_dict[loc][chan] = self[loc][chan].spectrogram(*args,
                                                                   **kwargs)
        return spec_dict

    def __copy__(self):
        new = type(self)()
        for key in self.keys():
            new[key] = self[key].copy()

    @classmethod
    def _gen_pwave(cls, stations, amplitude, phi, theta, frequency, duration,
                   Fs=100, c=3000, phase=0):
        """

        Parameters
        ----------
        stations : `seispy.station.StationArray` or `OrderedDict`
            ordered dict of stations, channels, locations
        amplitude : `float`
            amplitude of pwave
        phi
        theta
        frequency
        duration
        Fs
        c
        noise_amp
        phase
        segdur

        Returns
        -------

        """
        cphi = np.cos(phi)
        sphi = np.sin(phi)
        ctheta = np.cos(theta)
        stheta = np.sin(theta)
        src_dir = np.array([cphi * stheta, sphi * stheta, ctheta])
        # get time delays
        taus = np.array([-np.dot(src_dir, stations[key]) / c for key in
                         stations.keys()])
        tau_round = np.round(taus * Fs) / Fs
        ts = min(-tau_round)
        te = max(-tau_round)
        times = np.arange(0, duration, 1 / Fs)
        Nsamps = int(duration * Fs)
        # shift backward in time
        times += ts
        data = SeismometerArray()
        final_times = np.arange(0, duration, 1 / Fs)
        for key in stations.keys():
            data[key] = {}
            station = stations[key]
            delay = -np.dot(src_dir, station) / c
            # delaySamps = int(ts * Fs + np.round(delay * Fs))
            signal = np.zeros(times.size)
            if frequency == 0:
                signal = amplitude * np.random.randn(times.size)
            else:
                # most of our noise spectra will be one-sided, but this
                # is a real
                # signal, so we multiply this by two.
                signal = amplitude * np.sin(2 * np.pi * frequency *
                                            (times + delay) + phase)
            # impose time delay
            data[key]['HHE'] = Trace(src_dir[0] * signal, sample_rate=Fs,
                                     times=final_times, unit=u.m)
            data[key]['HHN'] = Trace(src_dir[1] * signal, sample_rate=Fs,
                                     times=final_times, unit=u.m)
            data[key]['HHZ'] = Trace(src_dir[2] * signal, sample_rate=Fs,
                                     times=final_times, unit=u.m)
            for key2 in data[key].keys():
                data[key][key2].location = station
        return data

    @classmethod
    def _gen_swave(cls, stations, amplitude, phi, theta, psi, frequency,
                   duration, phase=0, Fs=100, c=3000):
        """
        simulate s-wave in a certain direction

        Parameters
        ----------
        stations : `dict`
            dictionary of station locations
        amplitude : `float`
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
        src_dir = np.array([cphi * stheta, sphi * stheta, ctheta])
        # Get relative amplitudes in E,N,Z directions
        # based on polarizations. See internal method below.
        dx, dy, dz = get_polarization_coeffs(phi, theta, psi)

        # get time delays
        taus = np.array([-np.dot(src_dir, stations[key]) / c for key in
                         stations.keys()])
        tau_round = np.round(taus * Fs) / Fs
        ts = min(-tau_round)
        # te = max(-tau_round)
        # Nsamps = int(duration * Fs)
        final_times = np.arange(0, duration, 1 / Fs)
        times = np.arange(0, duration, 1 / Fs)
        # shift backward in time
        times += ts
        data = SeismometerArray()
        for key in stations.keys():
            data[key] = {}
            station = stations[key]
            delay = -np.dot(src_dir, station) / c
            signal = np.zeros(times.size)
            if frequency == 0:
                signal = amplitude * np.random.randn(times.size)
            else:
                signal = amplitude * np.sin(2 * np.pi * frequency *
                                            (times + delay) + phase)
            # impose time delay
            data[key]['HHE'] = Trace(dx * signal, sample_rate=Fs,
                                     times=final_times,
                                     unit=u.m, name=key)
            data[key]['HHN'] = Trace(dy * signal, sample_rate=Fs,
                                     times=final_times,
                                     unit=u.m, name=key)
            data[key]['HHZ'] = Trace(dz * signal, sample_rate=Fs,
                                     times=final_times,
                                     unit=u.m, name=key)
            for key2 in data[key].keys():
                data[key][key2].location = station
        return data

    @classmethod
    def _gen_rwave(cls, stations, amplitude, phi, theta, rwave_params,
                   frequency, duration, Fs=100, phase=0):
        """

        Parameters
        ----------
        stations
        amplitude
        phi
        theta
        epsilon
        alpha
        frequency
        duration
        Fs
        c
        noise_amp
        phase
        segdur

        Returns
        -------

        """
        # read r-wave eigenfunctions from dict
        C2 = rwave_params['C2']
        C4 = rwave_params['C4']
        Nvh = rwave_params['Nvh']
        a1 = rwave_params['a1']
        a2 = rwave_params['a2']
        a3 = rwave_params['a3']
        a4 = rwave_params['a4']
        v = rwave_params['v']

        # sines and cosines of source direction
        if theta != np.pi / 2:
            print('WARNING: Injecting R-Wave and theta!=pi/2')
        cphi = np.cos(phi)
        sphi = np.sin(phi)
        ctheta = np.cos(theta)
        stheta = np.sin(theta)

        # unit vector in direction of source direction
        src_dir = np.array([cphi * stheta, sphi * stheta, ctheta])

        # initialize data and times data structures
        data = SeismometerArray()
        times = np.arange(0, duration, 1 / Fs)
        for key in stations.keys():
            data[key] = {}
            station = stations[key]
            # time delay for this station
            delay = -np.dot(src_dir, station) / v
            # round time delay to nearest sampled time
            delay = np.round(delay * Fs) / Fs
            # initialize signal vector
            signal = np.zeros(times.size)
            if frequency == 0:
                signal = amplitude * np.random.randn(times.size)
            else:
                # most of our noise spectra will be one-sided,
                # but this is a real
                # signal, so we multiply this by two.
                signal = amplitude * np.cos(2 * np.pi * frequency *
                                            (times + delay) + phase)
                signal_phaseoff = amplitude * np.sin(2 * np.pi * frequency *
                                                     (times + delay) + phase)
            # compute eigenfunctions for this station
            k = 2 * np.pi * frequency / v
            z = station[2]

            r1 = ((np.exp(-a1 * k * z) + C2 * np.exp(-a2 * k * z)) *
                  (1 / (1 + C2)))
            r2 = ((np.exp(-a3 * k * z) + C4 * np.exp(-a4 * k * z)) *
                  (Nvh / (1 + C4)))

            # compute output of each channel
            data[key]['HHE'] = cphi * r1 * Trace(signal, sample_rate=Fs,
                                                 times=times, unit=u.m)
            data[key]['HHN'] = sphi * r1 * Trace(signal, sample_rate=Fs,
                                                 times=times, unit=u.m)
            data[key]['HHZ'] = -r2 * Trace(signal_phaseoff, sample_rate=Fs,
                                           times=times, unit=u.m)
            for key2 in data[key].keys():
                data[key][key2].location = station
        return data

    def add_p_wave(self, amplitude, phi, theta, frequency,
                   duration, phase=0, Fs=100, c=5700):
        """
        Add p-wave to seismometer array's data.
        Updates seismometer array data in place.

        Parameters
        ----------
        amplitude : `float`
            amplitude of p wave
        phi : `float`
            phi for p-wave injection
        theta : `float`
            polar angle for p-wave injection
        frequency : `float`
            frequency of injected wave
        duration : `float`
            duration of injected wave
        phase : `float`
            phase of injected sine-wave
        Fs : `float`, optional
            sample rate of wave we're generating
        c : `float`, optional
            velocity of injected wave
        """
        locations = self.get_locations()
        p_data = SeismometerArray._gen_pwave(locations, amplitude, phi, theta,
                                             frequency, duration,
                                             phase=phase, Fs=Fs, c=c
                                             )
        self.add_another_seismometer_array(p_data)

    def add_s_wave(self, amplitude, phi, theta, psi, frequency,
                   duration, phase=0, Fs=100, c=3000):
        """
        Add s-wave to seismometer array's data.
        Updates seismometer array data in place.

        Parameters
        ----------
        amplitude : `float`
            amplitude of shear wave
        phi : `float`
            phi for s-wave injection
        theta : `float`
            polar angle for s-wave injection
        psi : `float`
            polarization angle
        frequency : `float`
            frequency of injected wave
        duration : `float`
            duration of injected wave
        phase : `float`
            phase of injected sine-wave
        Fs : `float`, optional
            sample rate of wave we're generating
        c : `float`, optional
            velocity of injected wave
        """

        locations = self.get_locations()
        s_data = SeismometerArray._gen_swave(locations, amplitude, phi,
                                             theta, psi, frequency, duration,
                                             phase=phase, Fs=Fs, c=c
                                             )
        self.add_another_seismometer_array(s_data)

    def add_r_wave(self, amplitude, phi, theta, frequency,
                   duration, rayleigh_paramfile=None,
                   rayleigh_paramdict=None,
                   phase=0, Fs=100):
        """
        Add rayleigh wave to seismometer array's data.
        Updates seismometer array data in place.

        Parameters
        ----------
        amplitude : `float`
            amplitude of rayleigh wave
        phi : `float`
            phi for r-wave injection
        theta : `float`
            polar angle for r-wave injection
        epsilon : `float`
            vertical to horizontal amplitude ratio
        alpha : `float`
            depth attenuation factor
        frequency : `float`
            frequency of injected wave
        duration : `float`
            duration of injected wave
        phase : `float`
            phase of injected sine-wave
        Fs : `float`
            sample rate of wave we're generating
        c : `float`
            velocity of injected wave
        """
        if rayleigh_paramfile is None and rayleigh_paramdict is None:
            print('WARNING: No Rayleigh paramfile specified for injection')
            print('\tusing default eigenfunction')
            rwave_params = {'a1': 0.47,
                            'a3': 0.73,
                            'a2': 1.51,
                            'a4': 0.25,
                            'v': 2504,
                            'C2': -1.29,
                            'C4': 2.29,
                            'Nvh': -0.68}
        elif rayleigh_paramfile is not None:
            rwave_params = np.load(rayleigh_paramfile)[0]
        else:
            rwave_params = rayleigh_paramdict

        locations = self.get_locations()
        r_data = SeismometerArray._gen_rwave(locations, amplitude, phi,
                                             theta, rwave_params,
                                             frequency, duration, phase=phase,
                                             Fs=Fs)
        self.add_another_seismometer_array(r_data)

    @classmethod
    def initialize_all_good(cls, location_dict, duration, chans_type='useful',
                            start_time=0, Fs=None):
        data = cls()
        for name in location_dict.keys():
            data[name] = \
                Seismometer.initialize_all_good(duration,
                                                chans_type=chans_type,
                                                start_time=start_time,
                                                location=location_dict[name],
                                                name=name, Fs=Fs)
        return data

    @classmethod
    def _gen_white_gaussian_noise(cls, station_names, psd_amp, sample_rate,
                                  duration, segdur=None, seed=None):
        if segdur is None:
            segdur = duration
        data = SeismometerArray()
        psd = FrequencySeries(psd_amp * np.ones(100000), df=1. / segdur)
        psd[0] = 0
        for station in station_names:
            data[station] = Seismometer.initialize_all_good(duration=segdur)
            # set data channels to white noise
            data[station]['HHE'] = gaussian.noise_from_psd(duration,
                                                           sample_rate, psd,
                                                           seed=seed,
                                                           name=station,
                                                           unit=u.m)
            data[station]['HHN'] = gaussian.noise_from_psd(duration,
                                                           sample_rate, psd,
                                                           seed=seed,
                                                           name=station,
                                                           unit=u.m)
            data[station]['HHZ'] = gaussian.noise_from_psd(duration,
                                                           sample_rate, psd,
                                                           seed=seed,
                                                           name=station,
                                                           unit=u.m)
        return data

    def get_locations(self):
        location_dir = StationArray()
        for seismometer in self.keys():
            location_dir[seismometer] = \
                self[seismometer]['HHE'].location
        return location_dir

    def add_white_noise(self, psd_amp, segdur=None, seed=0):
        """
        Add some white noise to your seismometer array.
        Each channel with have a different noise realization
        but the same white noise amplitude.
        This done in place.
        """
        sensors = self.keys()
        # get duration
        Fs = self[sensors[0]]['HHE'].sample_rate.value
        duration = self[sensors[0]]['HHE'].size / Fs
        WN_array = \
            SeismometerArray._gen_white_gaussian_noise(sensors, psd_amp,
                                                       Fs, duration,
                                                       segdur=segdur,
                                                       seed=seed)
        self.add_another_seismometer_array(WN_array)

    def add_another_seismometer_array(self, other):
        """
        Internal method for adding two arrays
        that have all of the same station names.
        This shoudl only be used for simulation purposes.
        This is for combining noise and/or signal
        injections
        """
        for sensor in self.keys():
            # only care about HHE, HHZ, HHN for this
            # since we're simulating stuff
            self[sensor]['HHE'] += other[sensor]['HHE']
            self[sensor]['HHN'] += other[sensor]['HHN']
            self[sensor]['HHZ'] += other[sensor]['HHZ']

    def get_coherences_bband(self, flow, fhigh, channels=None, fftlength=2,
                             overlap=1, window='hann', nproc=8):
        """TODO: Docstring for get_coherences.

        Parameters
        ----------
        recovery_Freq : TODO

        Returns
        -------
        TODO

        """
        stations = self.keys()
        if channels is None:
            channels = ['HHE', 'HHN', 'HHZ']
        Ys = []
        ct = 0
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ii == jj and ll < kk:
                            # don't double count channels
                            continue
                        else:
                            # get csd spectrogram
                            self[station1][channels[kk]].bandpass(flow, fhigh)
                            self[station2][channels[ll]].bandpass(flow, fhigh)
                            P12 = \
                                self[station1][channels[kk]].csd_spectrogram(self[station2][channels[ll]],
                                                                             stride=fftlength,
                                                                             window=window,
                                                                             overlap=overlap,
                                                                             nproc=nproc)
                            P11 = \
                                self[station1][channels[kk]].csd_spectrogram(self[station1][channels[kk]],
                                                                             stride=fftlength,
                                                                             window=window,
                                                                             overlap=overlap,
                                                                             nproc=nproc)
                            P22 = \
                                self[station2][channels[ll]].csd_spectrogram(self[station2][channels[ll]],
                                                                             stride=fftlength,
                                                                             window=window,
                                                                             overlap=overlap,
                                                                             nproc=nproc)

                            P12 = P12.mean(0)
                            P11 = P11.mean(0)
                            P22 = P22.mean(0)
                            Y_of_f = P12/np.sqrt(P11 * P22)
                            Ys.append(FrequencySeries(Y_of_f, x0=0,
                                    dx=1./fftlength))
                            ct += 1
        return Ys


    def get_coherences(self, recovery_freq, channels=None,
                       stride=4, fftlength=2,
                       overlap=1, window='hann', nproc=8):
        """TODO: Docstring for get_coherences.

        Parameters
        ----------
        recovery_Freq : TODO

        Returns
        -------
        TODO

        """
        if stride < fftlength:
            stride = fftlength
        stations = self.keys()
        if channels is None:
            channels = ['HHE', 'HHN', 'HHZ']
        First = 1
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ii == jj and ll < kk:
                            # don't double count channels
                            continue
                        else:
                            # get csd spectrogram
                            P12 = \
                                self[station1][channels[kk]].csd_spectrogram(self[station2][channels[ll]],
                                                                             stride=stride,
                                                                             fftlength=fftlength,
                                                                             window=window,
                                                                             overlap=overlap,
                                                                             nproc=nproc)
                            # take mean across time
                            cp = P12.mean(0)
                            idx = \
                                np.where(cp.frequencies.value == float(recovery_freq))[0]
                            # multiply by 2
                            Y_of_t = (((P12[:, idx[0] - 2:idx[0] + 3]).sum(axis=1))).value * 2
                            if First:
                                Ys = Y_of_t.T
                                First = 0
                            else:
                                Ys = np.vstack((Ys, Y_of_t.T))
        return Ys

    def get_response_matrix_healpy(self, rec_type, station_locs, recovery_freq,
                                   v, nside=8, autocorrelations=True,
                                   paramfile=None, channels=None,
                                   fftlength=2, overlap=1,
                                   nproc=1, iter_lim=1000,
                                   atol=1e-6, btol=1e-6):
        """TODO: Docstring for get_gamma_matrices.
        Returns
        -------
        TODO

        """
        stations = self.keys()
        if channels is None:
            channels = ['HHE', 'HHN', 'HHZ']
        origin_vec = station_locs[stations[0]]
        First = True
        First = 1
        for ii, station1 in enumerate(stations):
            for kk, chan1 in enumerate(channels):
                phis = []
                thetas = []
                # get gamma matrix
                # convert freq to float
                rf = float(recovery_freq)
                # get diffraction limited spot size
                # get nside for healpy based on npix (taken up to
                # the next power of two)
                npix2 = hp.nside2npix(nside)
                # get theta and phi
                thetas_temp = np.zeros(npix2)
                phis_temp = np.zeros(npix2)
                for mm in range(npix2):
                    thetas_temp[mm], phis_temp[mm] =\
                        hp.pixelfunc.pix2ang(nside, mm)
                # get overlap reduction functions now
                if rec_type is 's':
                    # get s orf
                    g1, g2, g1_s, g2_s = \
                        response_picker(rec_type,
                                        set_channel_vector(channels[ll]),
                                        station_locs[station1],
                                        origin_vec, v,
                                        float(recovery_freq),
                                        thetas=thetas_temp,
                                        phis=phis_temp,
                                        epsilon=epsilon,
                                        alpha=alpha,
                                        healpy=True)
                    # append new, flattened, g onto the end
                    # of already generated on
                    g = np.vstack((g1, g2))
                    shapes = [g1_s, g2_s]
                    phis = np.hstack((phis_temp, phis_temp))
                    thetas = np.hstack((thetas_temp, thetas_temp))
                else:
                    # get p or r orf
                    g1, g_s = \
                        response_picker(rec_type, set_channel_vector(channels[kk]),
                                        station_locs[station1],
                                        origin_vec,
                                        v,
                                        float(recovery_freq),
                                        thetas=thetas_temp,
                                        phis=phis_temp,
                                        epsilon=epsilon,
                                        alpha=alpha,
                                        healpy=True)
                    # append new, flattened, g onto the
                    # end of the one we've generated
                    g = g1
                    phis = phis_temp
                    thetas = thetas_temp
                    shapes = [g_s]
                if First:
                    # for now for G:
                    # columns = channels
                    # rows = directions
                    G = g
                    First = 0
                else:
                    G = np.hstack((G, g))
        return G, phis, thetas, shapes

    def get_gamma_matrix_nonhealpy(self, rec_type, station_locs,
                                   v, recovery_freq, thetas, phis,
                                   autocorrelations=True,
                                   rayleigh_paramfile=None,
                                   rayleigh_paramdict=None,
                                   channels=None
                                   ):
        """TODO: Docstring for get_gamma_matrices.
        Returns
        -------
        TODO

        """
        stations = self.keys()
        if channels is None:
            channels = ['HHE', 'HHN', 'HHZ']
        # get theta and phi
        npairs = 0
        nstreams = (np.size(stations) * np.size(channels))
        npairs = (nstreams * (nstreams + 1)) / 2.
        if rec_type == 's':
            G = np.zeros((2 * thetas.size * phis.size, int(npairs)),
                         dtype=complex)
        else:
            G = np.zeros((thetas.size * phis.size, int(npairs)),
                         dtype=complex)

        ct = 0
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ii == jj and ll < kk:
                            # don't double count channels
                            continue
                        else:
                            # get gamma matrix
                            # convert freq to float
                            # get diffraction limited spot size
                            # get nside for healpy based on npix (taken up to
                            # the next power of two)
                            # get overlap reduction functions now
                            if rec_type is 's':
                                # get s orf
                                g1, g2, g1_s, g2_s = orf_picker(rec_type, set_channel_vector(channels[kk]),
                                                                set_channel_vector(channels[ll]),
                                                                station_locs[station1],
                                                                station_locs[station2], v,
                                                                float(recovery_freq),
                                                                thetas=thetas,
                                                                phis=phis,
                                                                healpy=False)
                                # append new, flattened, g onto the end
                                # of already generated on
                                g = np.vstack((g1, g2))
                                shapes = [g1_s, g2_s]
                            else:
                                # get p or r orf
                                g1, g_s = orf_picker(rec_type, set_channel_vector(channels[kk]),
                                                     set_channel_vector(channels[ll]), station_locs[station1],
                                                     station_locs[station2],
                                                     v,
                                                     float(recovery_freq),
                                                     thetas=thetas,
                                                     phis=phis,
                                                     rayleigh_paramfile=rayleigh_paramfile,
                                                     healpy=False)
                                # append new, flattened, g onto the
                                # end of the one we've generated
                                g = g1
                                #phis = phis_temp
                                #thetas = thetas_temp
                                shapes = [g_s]
                            G[:, ct] = g.squeeze()
                            ct += 1
        return G, shapes

    def get_gamma_matrix_healpy(self, rec_type, station_locs, v,
                                recovery_freq,
                                autocorrelations=True, rayleigh_paramfile=None,
                                rayleigh_paramdict=None,
                                channels=None, fftlength=2, overlap=1,
                                nproc=1, iter_lim=1000, atol=1e-6, btol=1e-6,
                                nside=8):
        """TODO: Docstring for get_gamma_matrices.
        Returns
        -------
        TODO

        """
        stations = self.keys()
        if channels is None:
            channels = ['HHE', 'HHN', 'HHZ']
        # get theta and phi
        npix2 = hp.nside2npix(nside)
        thetas_temp = np.zeros(npix2)
        phis_temp = np.zeros(npix2)
        pix = np.arange(npix2)

        thetas_temp, phis_temp =\
            hp.pixelfunc.pix2ang(nside, pix)
        npairs = 0
        nstreams = (np.size(stations) * np.size(channels))
        npairs = (nstreams * (nstreams + 1)) / 2.
        if rec_type == 's':
            G = np.zeros((2 * npix2, int(npairs)), dtype=complex)
        else:
            G = np.zeros((npix2, int(npairs)), dtype=complex)

        ct = 0
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ii == jj and ll < kk:
                            # don't double count channels
                            continue
                        else:
                            # get gamma matrix
                            # convert freq to float
                            # get diffraction limited spot size
                            # get nside for healpy based on npix (taken up to
                            # the next power of two)
                            # get overlap reduction functions now
                            if rec_type is 's':
                                # get s orf
                                g1, g2, g1_s, g2_s =\
                                    orf_picker(rec_type,
                                               set_channel_vector(channels[kk]),
                                               set_channel_vector(channels[ll]),
                                               station_locs[station1],
                                               station_locs[station2], v,
                                               float(recovery_freq),
                                               thetas=thetas_temp,
                                               phis=phis_temp,
                                               healpy=True)
                                # append new, flattened, g onto the end
                                # of already generated on
                                g = np.vstack((g1, g2))
                                shapes = [g1_s, g2_s]
                            else:
                                # get p or r orf
                                g1, g_s =\
                                    orf_picker(rec_type,
                                               set_channel_vector(channels[kk]),
                                               set_channel_vector(channels[ll]),
                                               station_locs[station1],
                                               station_locs[station2],
                                               v,
                                               float(recovery_freq),
                                               thetas=thetas_temp,
                                               phis=phis_temp,
                                               rayleigh_paramfile=rayleigh_paramfile,
                                               healpy=True)
                                # append new, flattened, g onto the
                                # end of the one we've generated
                                g = g1
                                shapes = [g_s]
                            G[:,ct] = g.squeeze()
                            ct += 1
        return G, shapes
