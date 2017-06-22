from __future__ import division
from collections import OrderedDict
from warnings import warn
from ..utils import *
import astropy.units as u
from ..noise import gaussian
from ..trace import Trace, fetch
from ..recoverymap import RecoveryMap
import numpy as np
from scipy.sparse.linalg import lsqr
from scipy.linalg import svd, pinv, pinv2, cholesky, lu, lu_factor, eigh
from gwpy.frequencyseries import FrequencySeries
from .station import StationArray, homestake
import geographiclib.geodesic as ggg
import healpy as hp


class Seismometer(OrderedDict):
    """
    Station data
    """

    @classmethod
    def fetch_data(cls, station_name, st, et, framedir='./', chans_type='useful'):
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
            seismometer[chan] = fetch(st, et, station_name + ':' + chan, framedir=framedir)
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
                                   sample_rate=Fs * u.Hz, epoch=start_time, name=name + ' East',
                                   unit=u.m, channel='%s:%s' % (name, 'HHE'))
        seismometer['HHN'] = Trace(np.zeros(int(duration * Fs)),
                                   sample_rate=Fs * u.Hz, epoch=start_time, name=name + ' North',
                                   unit=u.m, channel='%s:%s' % (name, 'HHN'))
        seismometer['HHZ'] = Trace(np.zeros(int(duration * Fs)),
                                   sample_rate=Fs * u.Hz, epoch=start_time, name=name + ' Vertical',
                                   unit=u.m, channel='%s:%s' % (name, 'HHZ'))
        if chans_type == 'fast_chans':
            for chan in seismometer.keys():
                seismometer[chan].location = location
            return seismometer
        else:
            seismometer['LCQ'] = Trace(100 * np.ones(int(duration * 1)),
                                       sample_rate=1 * u.Hz, epoch=start_time, name=name + ' Clock Quality',
                                       channel='%s:%s' % (name, 'LCQ'))
            seismometer['LCE'] = Trace(np.zeros(int(duration * 1)),
                                       sample_rate=1 * u.Hz, epoch=start_time, name=name + ' Clock Phase\
                    Error', channel='%s:%s' % (name, 'LCE'))
            seismometer['VM1'] = Trace(np.zeros(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz, epoch=start_time, name=name + ' Mass\
                    Position Channel 1', channel='%s:%s' % (name, 'VM1'))
            seismometer['VM2'] = Trace(np.zeros(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz, epoch=start_time, name=name + ' Mass\
                    Position Channel 2', channel='%s:%s' % (name, 'VM2'))
            seismometer['VM3'] = Trace(np.zeros(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz, epoch=start_time, name=name + ' Mass\
                    Position Channel 3', channel='%s:%s' % (name, 'VM3'))
            seismometer['VEP'] = Trace(13 * np.ones(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz, epoch=start_time, name=name + ' System\
                    Voltage', channel='%s:%s' % (name, 'VEP'))
            seismometer['VKI'] = Trace(np.zeros(int(duration * 0.1)),
                                       sample_rate=0.1 * u.Hz, epoch=start_time, name=name + ' System\
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
            seismometer with new channels for transverse and radial related to a single
            source
        """
        theta = np.radians(bearing)
        dataT = self['HHN'].value * np.cos(theta) - self['HHN'].value * np.sin(theta)
        dataR = self['HHE'].value * np.sin(theta) + self['HHE'].value * np.cos(theta)
        sample_rate = self['HHZ'].sample_rate
        epoch = self['HHZ'].epoch
        self['HHT'] = Trace(dataT, sample_rate=sample_rate, epoch=epoch, name='Transverse relative to %4.2f' % bearing,
                            channel=self['HHN'].channel)
        self['HHR'] = Trace(dataR, sample_rate=sample_rate, epoch=epoch, name='Radial relative to %4.2f' % bearing,
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
                                                  framedir=framedir, chans_type=chans_type)
        return arr

    @classmethod
    def _gen_pwave(cls, stations, amplitude, phi, theta, frequency, duration, Fs=100, c=3000,
                   phase=0):
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
        times = np.arange(0, np.abs(ts) + duration + te, 1 / Fs)
        Nsamps = int(duration * Fs)
        # shift backward in time
        times += ts
        data = SeismometerArray()
        final_times = np.arange(0, duration, 1 / Fs)
        for key in stations.keys():
            data[key] = {}
            station = stations[key]
            delay = -np.dot(src_dir, station) / c
            delaySamps = int(ts * Fs + np.round(delay * Fs))
            signal = np.zeros(times.size)
            if frequency == 0:
                signal = amplitude * np.random.randn(times.size)
            else:
                # most of our noise spectra will be one-sided, but this is a real
                # signal, so we multiply this by two.
                signal = amplitude * np.sin(2 * np.pi * frequency * times + phase)
            # impose time delay
            amp = np.roll(signal, delaySamps)[:Nsamps]
            data[key]['HHE'] = Trace(src_dir[0] * amp, sample_rate=Fs,
                                     times=final_times, unit=u.m)
            data[key]['HHN'] = Trace(src_dir[1] * amp, sample_rate=Fs,
                                     times=final_times, unit=u.m)
            data[key]['HHZ'] = Trace(src_dir[2] * amp, sample_rate=Fs,
                                     times=final_times, unit=u.m)
            for key2 in data[key].keys():
                data[key][key2].location = station
        return data

    @classmethod
    def _gen_swave(cls, stations, amplitude, phi, theta, psi, frequency, duration,
                   phase=0, Fs=100, c=3000):
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
        te = max(-tau_round)
        Nsamps = int(duration * Fs)
        final_times = np.arange(0, duration, 1 / Fs)
        times = np.arange(0, np.abs(ts) + duration + te, 1 / Fs)
        # shift backward in time
        times += ts
        data = SeismometerArray()
        ct = 0
        for key in stations.keys():
            data[key] = {}
            station = stations[key]
            delay = -np.dot(src_dir, station) / c
            delaySamps = int(ts * Fs + np.round(delay * Fs))
            signal = np.zeros(times.size)
            if frequency == 0:
                signal = amplitude * np.random.randn(times.size)
            else:
                signal = amplitude * np.sin(2 * np.pi * frequency * times + phase)
            # impose time delay
            amp = np.roll(signal, delaySamps)[:Nsamps]
            data[key]['HHE'] = Trace(dx * amp, sample_rate=Fs, times=final_times,
                                     unit=u.m, name=key)
            data[key]['HHN'] = Trace(dy * amp, sample_rate=Fs, times=final_times,
                                     unit=u.m, name=key)
            data[key]['HHZ'] = Trace(dz * amp, sample_rate=Fs, times=final_times,
                                     unit=u.m, name=key)
            for key2 in data[key].keys():
                data[key][key2].location = station
        return data

    @classmethod
    def _gen_rwave(cls, stations, amplitude, phi, theta, epsilon,
                   alpha, frequency, duration, Fs=100, c=3000,
                   phase=0):
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
        times = np.arange(0, np.abs(ts) + duration + te, 1 / Fs)
        Nsamps = int(duration * Fs)
        # shift backward in time
        times += ts
        data = SeismometerArray()
        ct = 0
        final_times = np.arange(0, duration, 1 / Fs)
        for key in stations.keys():
            data[key] = {}
            station = stations[key]
            delay = -np.dot(src_dir, station) / c
            delaySamps = int(ts * Fs + np.round(delay * Fs))
            signal = np.zeros(times.size)
            if frequency == 0:
                signal = amplitude * np.random.randn(times.size)
            else:
                # most of our noise spectra will be one-sided, but this is a real
                # signal, so we multiply this by two.
                signal = amplitude * np.cos(2 * np.pi * frequency * times + phase)
                signal_phaseoff = -amplitude * np.sin(2 * np.pi * frequency * times + phase)
            # impose time delay
            amp = np.roll(signal, delaySamps)[:Nsamps] * np.exp(-station[2] / alpha)
            amp2 = np.roll(signal_phaseoff, delaySamps)[:Nsamps] * \
                   np.exp(station[2] / alpha)
            data[key]['HHE'] = cphi * Trace(amp, sample_rate=Fs,
                                            times=final_times, unit=u.m)
            data[key]['HHN'] = sphi * Trace(amp, sample_rate=Fs,
                                            times=final_times, unit=u.m)
            data[key]['HHZ'] = epsilon * Trace(amp2, sample_rate=Fs,
                                               times=final_times, unit=u.m)
            for key2 in data[key].keys():
                data[key][key2].location = station
        return data

    def add_p_wave(self, amplitude, phi, theta, frequency,
                   duration, phase=0, Fs=100, c=5700):
        """
        Add s-wave to seismometer array's data. Updates seismometer array data in place.

        Parameters
        ----------
        amplitude : `float`
            amplitude of rayleigh wave
        phi : `float`
            phi for r-wave injection
        theta : `float`
            polar angle for r-wave injection
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
                                             frequency, duration, phase=phase, Fs=Fs, c=c
                                             )
        self.add_another_seismometer_array(p_data)

    def add_s_wave(self, amplitude, phi, theta, psi, frequency,
                   duration, phase=0, Fs=100, c=3000):
        """
        Add s-wave to seismometer array's data. Updates seismometer array data in place.

        Parameters
        ----------
        amplitude : `float`
            amplitude of rayleigh wave
        phi : `float`
            phi for r-wave injection
        theta : `float`
            polar angle for r-wave injection
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
                                             theta, psi, frequency, duration, phase=phase, Fs=Fs, c=c
                                             )
        self.add_another_seismometer_array(s_data)

    def add_r_wave(self, amplitude, phi, theta, epsilon, alpha, frequency,
                   duration, phase=0, Fs=100, c=200):
        """

        Add rayleigh wave to seismometer array's data. Updates seismometer array data in place.

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
        locations = self.get_locations()
        r_data = SeismometerArray._gen_rwave(locations, amplitude, phi,
                                             theta, epsilon, alpha, frequency, duration, phase=phase, Fs=Fs, c=c
                                             )
        self.add_another_seismometer_array(r_data)

    @classmethod
    def initialize_all_good(cls, location_dict, duration, chans_type='useful',
                            start_time=0, Fs=None):
        data = cls()
        for name in location_dict.keys():
            data[name] = Seismometer.initialize_all_good(duration,
                                                         chans_type=chans_type, start_time=start_time,
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
                                                           sample_rate, psd, seed=seed, name=station, unit=u.m)
            data[station]['HHN'] = gaussian.noise_from_psd(duration,
                                                           sample_rate, psd, seed=seed, name=station, unit=u.m)
            data[station]['HHZ'] = gaussian.noise_from_psd(duration,
                                                           sample_rate, psd, seed=seed, name=station, unit=u.m)
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
        WN_array = SeismometerArray._gen_white_gaussian_noise(sensors, psd_amp,
                                                              Fs, duration, segdur=segdur, seed=seed)
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

    def recovery_matrices(self, rec_str, station_locs, recovery_freq,
                          v_list, autocorrelations=True, epsilon=0.1, alpha=1000,
                          channels=None, phis=None, thetas=None, fftlength=2, overlap=1,
                          nproc=1, iter_lim=1000, atol=1e-6, btol=1e-6):
        """
        Recovery calculated with gamma*Y, gammaT*gamma into lsqr (the typical
        way we do it)

        Parameters
        ----------
        rec_str
        station_locs
        recovery_freq
        v_list
        autocorrelations
        epsilon
        alpha
        channels
        phis
        thetas
        fftlength
        overlap
        nproc
        iter_lim
        atol
        btol

        Returns
        -------

        """
        stations = self.keys()
        if channels is None:
            channels = ['HHE', 'HHN', 'HHZ']
        First = True
        distances = []
        for station in stations:
            for station2 in stations:
                diff_dir = station_locs[station] - station_locs[stations2]
                distances.append(np.sqrt(np.dot(diff_dir,diff_dir)))
        maxd = np.max(distances)

        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ll < kk:
                            # don't double count channels
                            continue
                        else:
                            # get csd spectrogram
                            P12 = \
                                self[station1][channels[kk]].csd_spectrogram(self[station2][channels[ll]],
                                                                             stride=fftlength,
                                                                             window='hann', overlap=overlap,
                                                                             nproc=nproc)
                            # take mean across time
                            cp = P12.mean(0)
                            idx = \
                                np.where(cp.frequencies.value == float(recovery_freq))[0]
                            # add up the frequency bin we want and adjacent ones
                            # to account for spectral leakage
                            p12 = cp[idx[0] - 1:idx[0] + 2].sum()
                            # reset g
                            g = []
                            shapes = []
                            for rec, v in zip(rec_str, v_list):
                                spot_size = v / (2 * maxd * recovery_freq)
                                npix = 4*np.pi / (spot_size **2)
                                if rec is 's':
                                    # get s orf
                                    g1, g2, g1_s, g2_s = orf_picker(rec, set_channel_vector(channels[kk]),
                                                                    set_channel_vector(channels[ll]),
                                                                    station_locs[station1],
                                                                    station_locs[station2], v,
                                                                    float(recovery_freq), thetas=thetas, phis=phis,
                                                                    epsilon=epsilon, alpha=alpha)
                                    shapes.append(g1_s)
                                    shapes.append(g2_s)
                                    # append new, flattened, g onto the end
                                    # of already generated on
                                    if len(g) > 0:
                                        g = np.vstack((g, g1, g2))
                                    else:
                                        g = np.vstack((g1, g2))
                                else:
                                    # get p or r orf
                                    g1, g_s = orf_picker(rec, set_channel_vector(channels[kk]),
                                                         set_channel_vector(channels[ll]), station_locs[station1],
                                                         station_locs[station2],
                                                         v, float(recovery_freq), thetas=thetas, phis=phis,
                                                         epsilon=epsilon, alpha=alpha)
                                    # append new, flattened, g onto the
                                    # end of the one we've generated
                                    try:
                                        g = np.vstack((g, g1))
                                    except ValueError:
                                        g = g1
                                    shapes.append(g_s)
                            # generate or add to gammaT*gamma or gammaT*Y
                            # vector
                            if First:
                                GG = np.dot(np.conj(g), np.transpose(g))
                                GY = np.conj(g) * p12
                                First = 0
                            else:
                                GG += np.dot(np.conj(g), np.transpose(g))
                                GY += np.conj(g) * p12
        # run least squares solver
        print 'IS GG SYMMETRIC??'
        print (GG.transpose() == GG).all()
        S = lsqr(np.real(GG), np.real(GY.value), iter_lim=iter_lim, atol=atol,
                 btol=btol, calc_var=True)
        R = S[-1] * np.real(GG)
        print R
        G_inv = np.linalg.pinv(np.real(GG))
        plt.pcolormesh(R, cmap='viridis')
        plt.colorbar()
        plt.savefig('resolution_matrix')
        maps = {}
        idx_low = 0
        # separate result into proper maps
        # with proper class
        if thetas is None:
            thetas = np.arange(3, 180, 6) * np.pi / 180
        if phis is None:
            phis = np.arange(3, 360, 6) * np.pi / 180
        for ii, rec in enumerate(rec_str):
            if rec is 's':
                length = shapes[ii][0] * shapes[ii][1]
                maps['s1'] = \
                    RecoveryMap(S[0].reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                                thetas, phis, 's1')
                maps['sigma_sq_s1'] =\
                    RecoveryMap(sig2.reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas,
                            phis, '\sigma^2_{s2}')
                idx_low += length
                maps['s2'] = \
                    RecoveryMap(S[0].reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                                thetas, phis, 's2')
                maps['sigma_sq_s2'] =\
                    RecoveryMap(S[-1].reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas,
                            phis, '\sigma^2_{s2}')
            else:
                length = shapes[ii][0] * shapes[ii][1]
                maps[rec] = \
                    RecoveryMap(S[0].reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                                thetas, phis, rec)
                maps['sigma_sq_'+rec] =\
                    RecoveryMap(S[-1].reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas,
                            phis, '\sigma^2_{'+rec+'}')
                idx_low += length

        print('Stopped at iteration number ' + str(S[2]))
        if S[1] == 1:
            print("We've found an exact solution")
        if S[1] == 2:
            print("We found an approximate solution")
        print('Converged to a relative residual of ' + str(S[3] /
                                                           np.sqrt((np.abs(GY.value) ** 2).sum())))
        return maps, phis, thetas
    def get_ffts(self, recovery_freq, channels=None, fftlength=2,
            overlap=1, window='hann',nproc=8):
        """TODO: Docstring for get_ffts.

        Parameters
        ----------
        recovery_Freq : TODO

        Returns
        -------
        TODO

        """
        stations = self.keys()
        if channels is None:
            channels = ['HHE','HHN','HHZ']
        First = 1
        data_samples = []
        for ii, station1 in enumerate(stations):
            for kk, channel in enumerate(channels):
                # get csd spectrogram
                fftgram_temp = \
                    fftgram(self[station1][channels[kk]],stride=fftlength
                                                         )
                # get covariance matrix entry
                cp = fftgram_temp.mean(0)
                idx = \
                    np.where(cp.frequencies.value == float(recovery_freq))[0]

                amax =  np.argmax(fftgram_temp[0,:])
                fft_of_t = fftgram_temp[:,idx[0]-1:idx[0]+2].sum(axis=1).value
                if First:
                    data_samples = fft_of_t.T
                    First = 0
                else:
                    data_samples = np.vstack((data_samples, fft_of_t.T))
        return data_samples


    def get_coherences(self, recovery_freq, channels=None, fftlength=2,
            overlap=1, window='hann',nproc=8):
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
            channels = ['HHE','HHN','HHZ']
        First = 1
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ii==jj and ll < kk:
                            # don't double count channels
                            continue
                        else:
                            # get csd spectrogram
                            P12 = \
                                self[station1][channels[kk]].csd_spectrogram(self[station2][channels[ll]],
                                                                             stride=fftlength,
                                                                             window=window,
                                                                             overlap=overlap,
                                                                             nproc=nproc)
                            # get covariance matrix entry
                            n = np.mean(P12*np.conj(P12), 0) -\
                                np.mean(P12,0)*np.mean(np.conj(P12), 0)

                            # take mean across time
                            cp = P12.mean(0)
                            idx = \
                                np.where(cp.frequencies.value == float(recovery_freq))[0]

                            Y_of_t = P12[:,idx[0]-1:idx[0]+2].sum(axis=1).value
                            if First:
                                Ys = Y_of_t.T
                                First = 0
                            else:
                                Ys = np.vstack((Ys, Y_of_t.T))
        return Ys

    def get_response_matrix_healpy(self, rec_type, station_locs, recovery_freq,
                          v, nside=8, autocorrelations=True, epsilon=0.1, alpha=1000,
                          channels=None, fftlength=2, overlap=1,
                          nproc=1, iter_lim=1000, atol=1e-6, btol=1e-6):
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
                        hp.pixelfunc.pix2ang(nside,mm)
                # get overlap reduction functions now
                if rec_type is 's':
                    # get s orf
                    g1, g2, g1_s, g2_s = response_picker(rec_type,
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
                    g1, g_s = response_picker(rec_type, set_channel_vector(channels[kk]),
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

    def get_gamma_matrix_healpy(self, rec_type, station_locs, recovery_freq,
                          v, autocorrelations=True, epsilon=0.1, alpha=1000,
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
        First = True
        distances = []
        for station in stations:
            for station2 in stations:
                diff_dir = np.asarray(station_locs[station]) -\
                    np.asarray(station_locs[station2])
                distances.append(np.sqrt(np.dot(diff_dir,diff_dir)))
        maxd = np.max(distances)
        First = 1
        npix2 = hp.nside2npix(nside)
        # get theta and phi
        thetas_temp = np.zeros(npix2)
        phis_temp = np.zeros(npix2)
        pix = np.arange(npix2)

        thetas_temp, phis_temp =\
        hp.pixelfunc.pix2ang(nside,pix)
        npairs = 0
        # get number of pairs of stations
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ii==jj and ll < kk:
                            # don't double count channels
                            continue
                        npairs += 1
#        G = np.zeros((npix2, npairs))
        ct = 0
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ll < kk:
                            # don't double count channels
                            continue
                        if ii == jj and ll == kk:
                            continue
                        else:
                            phis = []
                            thetas = []
                            # get gamma matrix
                            # convert freq to float
                            rf = float(recovery_freq)
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
                                g1, g_s = orf_picker(rec_type, set_channel_vector(channels[kk]),
                                                     set_channel_vector(channels[ll]), station_locs[station1],
                                                     station_locs[station2],
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
                                G = g
                                First = 0
                            else:
                                G = np.hstack((G,g))
                            ct += 1
        return G, phis, thetas, shapes

    def recovery_matrices_pinv(self, rec_str, station_locs, recovery_freq,
                          v_list, autocorrelations=True, epsilon=0.1, alpha=1000,
                          channels=None, phis=None, thetas=None, fftlength=2, overlap=1,
                          nproc=1, iter_lim=1000, atol=1e-6, btol=1e-6):
        """
        Calculate model parameters, S, using pseudo-inverse of gamma matrix

        Parameters
        ----------
        rec_str
        station_locs
        recovery_freq
        v_list
        autocorrelations
        epsilon
        alpha
        channels
        phis
        thetas
        fftlength
        overlap
        nproc
        iter_lim
        atol
        btol

        Returns
        -------

        """
        stations = self.keys()
        if channels is None:
            channels = ['HHE', 'HHN', 'HHZ']
        First = True
        ct = 0
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ll < kk:
                            # don't double count channels
                            continue
                        if ii == jj and ll == kk:
                            continue
                        else:
                            # get csd spectrogram
                            ct += 1
                            P12 = \
                                self[station1][channels[kk]].csd_spectrogram(self[station2][channels[ll]],
                                                                             stride=fftlength,
                                                                             window='hann',
                                                                             overlap=0,
                                                                             nproc=nproc)
                            # get covariance matrix entry
                            n = np.mean(P12*np.conj(P12), 0) -\
                                np.mean(P12,0)*np.mean(np.conj(P12), 0)

                            # take mean across time
                            cp = P12.mean(0)
                            idx = \
                                np.where(cp.frequencies.value == float(recovery_freq))[0]

                            Y_of_t = P12[:,idx[0]-1:idx[0]+2].sum(axis=1).value
                            # add up the frequency bin we want and adjacent ones
                            # to account for spectral leakage
                            p12 = cp[idx[0] - 1:idx[0] + 2].sum()
                            n = n[idx[0]-1:idx[0]+2].sum()
                            # reset g
                            g = []
                            shapes = []
                            for rec, v in zip(rec_str, v_list):
                                if rec is 's':
                                    # get s orf
                                    g1, g2, g1_s, g2_s = orf_picker(rec, set_channel_vector(channels[kk]),
                                                                    set_channel_vector(channels[ll]),
                                                                    station_locs[station1],
                                                                    station_locs[station2], v,
                                                                    float(recovery_freq), thetas=thetas, phis=phis,
                                                                    epsilon=epsilon, alpha=alpha)
                                    shapes.append(g1_s)
                                    shapes.append(g2_s)
                                    # append new, flattened, g onto the end
                                    # of already generated on
                                    if len(g) > 0:
                                        g = np.vstack((g, g1, g2))
                                    else:
                                        g = np.vstack((g1, g2))
                                else:
                                    # get p or r orf
                                    g1, g_s = orf_picker(rec, set_channel_vector(channels[kk]),
                                                         set_channel_vector(channels[ll]), station_locs[station1],
                                                         station_locs[station2],
                                                         v, float(recovery_freq), thetas=thetas, phis=phis,
                                                         epsilon=epsilon, alpha=alpha)
                                    # append new, flattened, g onto the
                                    # end of the one we've generated
                                    try:
                                        g = np.vstack((g, g1))
                                    except ValueError:
                                        g = g1
                                    shapes.append(g_s)
                            # generate or add to gammaT*gamma or gammaT*Y
                            # vector
                            if First:
                                # for now for G:
                                # columns = channels
                                # rows = directions
                                G = g
                                Ys = Y_of_t.T
                                N = [n.value]
                                First = 0
                            else:
                                G = np.hstack((G, g))
                                Ys = np.vstack((Ys, Y_of_t.T))
                                N.append(n.value)
        print Ys.shape
        print ct
        G = G.T
        # estimate of power in given direction
        Y = np.mean(Ys,1)
        X = np.real(np.dot(G.T, Y))
        N_inv = np.diag(np.ones(X.size))
        N = np.diag(np.ones(X.size))
        # covariance of estimate of power in each direction
        # this is a covariance across directions
        # ok now let's get an inverse of our covariance matrix
        # explicit
        # beam pattern function
        BigGamma = np.real(np.dot(G.T, G))
        inv_BigGamma = pinv2(BigGamma, rcond=1e-5)
        # Now let's use this to get estimate of power in each direction
        # where we're effectively deconvolving sky directions
        # This is pseudo-inverse of beam pattern matrix
        CleanMap = np.real(np.dot(inv_BigGamma, X))
        # our map is the clean map
        S = [CleanMap]
        # the variance is roughly the fisher matrix
        # times its inverse
        sig2 = np.real(np.dot(BigGamma, np.dot(N_inv, BigGamma.T.conj())))
        std_estimator = np.sqrt(sig2)
        s = X - np.dot(BigGamma, S[0])
        logL_n = -0.5 * np.dot(np.dot(X.conj().T, N_inv), X)
        plt.figure()
        plt.pcolormesh(N)
        plt.colorbar()
        plt.title('covariance')
        plt.savefig('covariance')
        plt.close()

        plt.figure()
        plt.pcolormesh(sig2)
        plt.title('sigma square')
        plt.colorbar()
        plt.savefig('sigma_sq')
        plt.close()

        plt.figure()
        plt.pcolormesh(BigGamma)
        plt.colorbar()
        plt.title('beam pattern')
        plt.savefig('beam_pattern')
        plt.close()
        # separate result into proper maps
        # with proper class
        idx_low = 0
        maps = {}
        if thetas is None:
            thetas = np.arange(3, 180, 6) * np.pi / 180
        if phis is None:
            phis = np.arange(3, 360, 6) * np.pi / 180
        for ii, rec in enumerate(rec_str):
            if rec is 's':
                length = shapes[ii][0] * shapes[ii][1]
                maps['s1'] = \
                    RecoveryMap(np.real(S[0]).reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                                thetas, phis, 's1')
                maps['sigma_sq_s1'] =\
                    RecoveryMap(sig2.reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas,
                            phis, '\sigma^2_{s2}')
                idx_low += length
                maps['s2'] = \
                    RecoveryMap(np.real(S[0]).reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                                thetas, phis, 's2')
                maps['sigma_sq_s2'] =\
                    RecoveryMap(np.real(np.diag(std_estimator)).reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas,
                            phis, '\sigma^2_{s2}')
            else:
                length = shapes[ii][0] * shapes[ii][1]
                maps[rec] = \
                    RecoveryMap(S[0].reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                                thetas, phis, rec)
                maps['sigma_sq_'+rec] =\
                    RecoveryMap(np.real(np.diag(std_estimator)).reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas,
                            phis, '\sigma_{'+rec+'}')
                maps['SNR'] =\
                    RecoveryMap(np.real(S[0]).reshape(g.shape)[idx_low:idx_low
                        + length].reshape(shapes[ii])/np.real(np.diag(std_estimator)).reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas,
                            phis, 'SNR')


                idx_low += length
        return maps, phis, thetas
    def recovery_matrices_pinv_healpy(self, rec_str, station_locs, recovery_freq,
                          v_list, autocorrelations=True, epsilon=0.1, alpha=1000,
                          channels=None, fftlength=2, overlap=1,
                          nproc=1, iter_lim=1000, atol=1e-6, btol=1e-6):
        """
        Calculate model parameters, S, using pseudo-inverse of gamma matrix

        Parameters
        ----------
        rec_str
        station_locs
        recovery_freq
        v_list
        autocorrelations
        epsilon
        alpha
        channels
        phis
        thetas
        fftlength
        overlap
        nproc
        iter_lim
        atol
        btol

        Returns
        -------

        """
        stations = self.keys()
        if channels is None:
            channels = ['HHE', 'HHN', 'HHZ']
        print channels
        First = True
        distances = []
        for station in stations:
            for station2 in stations:
                diff_dir = np.asarray(station_locs[station]) -\
                    np.asarray(station_locs[station2])
                distances.append(np.sqrt(np.dot(diff_dir,diff_dir)))
        maxd = np.max(distances)
        ct = 0
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ll < kk:
                            # don't double count channels
                            continue
                        if ii == jj and ll == kk:
                            continue
                        else:
                            # get csd spectrogram
                            ct +=1
                            P12 = \
                                self[station1][channels[kk]].csd_spectrogram(self[station2][channels[ll]],
                                                                             stride=fftlength,
                                                                             window='hann',
                                                                             overlap=0,
                                                                             nproc=nproc)
                            # get covariance matrix entry
                            n = np.mean(P12*np.conj(P12), 0) -\
                                np.mean(P12,0)*np.mean(np.conj(P12), 0)

                            # take mean across time
                            cp = P12.mean(0)
                            idx = \
                                np.where(cp.frequencies.value == float(recovery_freq))[0]

                            Y_of_t = P12[:,idx[0]-1:idx[0]+2].sum(axis=1).value
                            # add up the frequency bin we want and adjacent ones
                            # to account for spectral leakage
                            p12 = cp[idx[0] - 1:idx[0] + 2].sum()
                            n = n[idx[0]-1:idx[0]+2].sum()
                            # reset g
                            g = []
                            shapes = []
                            ang_shapes = []
                            thetas = []
                            phis = []
                            for rec, v in zip(rec_str, v_list):
                                rf = float(recovery_freq)
                                spot_size = v/(2*maxd*rf)
                                npix = (4 * np.pi / (spot_size **
                                    2))*np.sqrt(len(self.keys()))
                                nside = int(2 ** round(np.log2(round((int(npix)+1) **
                                    (0.5)))))
                                nside=4
                                npix2 = hp.nside2npix(nside)
                                thetas_temp = np.zeros(npix2)
                                phis_temp = np.zeros(npix2)
                                # get angles
                                for mm in range(npix2):
                                    thetas_temp[mm], phis_temp[mm] =\
                                        hp.pixelfunc.pix2ang(nside,mm)
                                if rec is 's':
                                    # get s orf
                                    g1, g2, g1_s, g2_s = orf_picker(rec, set_channel_vector(channels[kk]),
                                                                    set_channel_vector(channels[ll]),
                                                                    station_locs[station1],
                                                                    station_locs[station2], v,
                                                                    float(recovery_freq),
                                                                    thetas=thetas_temp,
                                                                    phis=phis_temp,
                                                                    epsilon=epsilon,
                                                                    alpha=alpha,
                                                                    healpy=True)
                                    shapes.append(g1_s)
                                    shapes.append(g2_s)
                                    # append new, flattened, g onto the end
                                    # of already generated on
                                    if len(g) > 0:
                                        g = np.vstack((g, g1, g2))
                                    else:
                                        g = np.vstack((g1, g2))
                                    ang_shapes.append(len(phis_temp))
                                else:
                                    # get p or r orf
                                    g1, g_s = orf_picker(rec, set_channel_vector(channels[kk]),
                                                         set_channel_vector(channels[ll]), station_locs[station1],
                                                         station_locs[station2],
                                                         v,
                                                         float(recovery_freq),
                                                         thetas=thetas_temp,
                                                         phis=phis_temp,
                                                         epsilon=epsilon,
                                                         alpha=alpha,
                                                         healpy=True)
                                    # append new, flattened, g onto the
                                    # end of the one we've generated
                                    try:
                                        g = np.vstack((g, g1))
                                    except ValueError:
                                        g = g1
                                    shapes.append(g_s)
                                    ang_shapes.append(len(phis_temp))
                            # generate or add to gammaT*gamma or gammaT*Y
                            # vector
                            if First:
                                # for now for G:
                                # columns = channels
                                # rows = directions
                                thetas = thetas_temp
                                phis = phis_temp
                                G = g
                                Ys = Y_of_t.T
                                N = [n.value]
                                First = 0
                            else:
                                phis = np.hstack((phis, phis_temp))
                                thetas = np.hstack((thetas, thetas_temp))
                                G = np.hstack((G, g))
                                Ys = np.vstack((Ys, Y_of_t.T))
                                N.append(n.value)
                            print phis.size
                            print phis_temp.size
        print ct
        print Ys.shape
###########################
###########################
###########################
###########################
        G = G.T
        # estimate of power in given direction
        Y = np.mean(Ys,1)
        X = np.real(np.dot(G.T, Y))
        N_inv = np.diag(np.ones(X.size))
        N = np.diag(np.ones(X.size))
        # covariance of estimate of power in each direction
        # this is a covariance across directions
        # ok now let's get an inverse of our covariance matrix
        # explicit
        # beam pattern function
        BigGamma = np.real(np.dot(G.T, G))
        inv_BigGamma = pinv2(BigGamma, rcond=1e-5)
        # Now let's use this to get estimate of power in each direction
        # where we're effectively deconvolving sky directions
        # This is pseudo-inverse of beam pattern matrix
        CleanMap = np.real(np.dot(inv_BigGamma, X))
        # our map is the clean map
        S = [CleanMap]
        # the variance is roughly the fisher matrix
        # times its inverse
        sig2 = np.real(np.dot(BigGamma, np.dot(N_inv, BigGamma.T.conj())))
        std_estimator = np.sqrt(sig2)
        s = X - np.dot(BigGamma, S[0])
        logL_n = -0.5 * np.dot(np.dot(X.conj().T, N_inv), X)
###########################
###########################
###########################
###########################
        idx_low = 0
        shape_idx_low = 0
        maps = {}
        print 'STOP'
        print phis.shape
        print thetas.shape
        print ang_shapes
        for ii, rec in enumerate(rec_str):
            phis_new = phis[shape_idx_low:shape_idx_low + ang_shapes[ii]]
            thetas_new = thetas[shape_idx_low:shape_idx_low + ang_shapes[ii]]
            if rec is 's':
                length = shapes[ii][0]
                maps['s1'] = \
                    RecoveryMap(np.real(S[0]).reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                            thetas_new, phis_new, 's1')
                maps['sigma_sq_s1'] =\
                    RecoveryMap(sig2.reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas_new,
                            phis_new, '\sigma^2_{s2}')
                idx_low += length
                maps['s2'] = \
                    RecoveryMap(np.real(S[0]).reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                            thetas_new, phis_new, 's2')
                maps['sigma_sq_s2'] =\
                    RecoveryMap(np.real(np.diag(std_estimator)).reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas_new,
                            phis_new, '\sigma^2_{s2}')
            else:
                length = shapes[ii][0]
                maps[rec] = \
                    RecoveryMap(S[0].reshape(g.shape)[idx_low:idx_low + length].reshape(shapes[ii]),
                            thetas_new, phis_new, rec)
                maps['sigma_sq_'+rec] =\
                        RecoveryMap(np.real(np.diag(std_estimator)).reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),
                                thetas_new,
                                phis_new, '\sigma_{'+rec+'}')
                maps['SNR'] =\
                    RecoveryMap(np.real(S[0]).reshape(g.shape)[idx_low:idx_low
                        +
                        length].reshape(shapes[ii])/np.real(np.diag(std_estimator)).reshape(g.shape)[idx_low:idx_low+length].reshape(shapes[ii]),thetas_new,
                        phis_new, 'SNR')

                idx_low += length
                shape_idx_low += ang_shapes[ii]
        return maps, phis, thetas

   def get_gamma_matrix_healpy(self, rec_type, station_locs, recovery_freq,
                          v, autocorrelations=True, epsilon=0.1, alpha=1000,
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
        First = True
        distances = []
        for station in stations:
            for station2 in stations:
                diff_dir = np.asarray(station_locs[station]) -\
                    np.asarray(station_locs[station2])
                distances.append(np.sqrt(np.dot(diff_dir,diff_dir)))
        maxd = np.max(distances)
        First = 1
        npix2 = hp.nside2npix(nside)
        # get theta and phi
        thetas_temp = np.zeros(npix2)
        phis_temp = np.zeros(npix2)
        pix = np.arange(npix2)

        thetas_temp, phis_temp =\
        hp.pixelfunc.pix2ang(nside,pix)
        npairs = 0 
        # get number of pairs of stations (AAM: this should be m*(m+1)/2 where m=Nstations*Nchannels )
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        if ii==jj and ll < kk:
                            # don't double count channels
                            continue
                        npairs += 1
#        G = np.zeros((npix2, npairs))
        ct = 0
        for ii, station1 in enumerate(stations):
            for jj, station2 in enumerate(stations):
                for kk, chan1 in enumerate(channels):
                    for ll, chan2 in enumerate(channels):
                        if jj < ii:
                            # we don't double count stations
                            continue
                        # AAM: commented out below (bug)
                        #if ll < kk:
                        #    # don't double count channels
                        #    continue
                        #if ii == jj and ll == kk:
                        if ii == jj and ll < kk:
                            continue
                        else:
                            phis = []
                            thetas = []
                            # get gamma matrix
                            # convert freq to float
                            rf = float(recovery_freq)
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
                                g1, g_s = orf_picker(rec_type, set_channel_vector(channels[kk]),
                                                     set_channel_vector(channels[ll]), station_locs[station1],
                                                     station_locs[station2],
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
                                G = g
                                First = 0
                            else:
                                G = np.hstack((G,g))
                            ct += 1
        return G, phis, thetas, shapes

