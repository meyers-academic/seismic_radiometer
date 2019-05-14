#! /usr/bin/python
import seispy.utils.responses as rsps
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from seispy import plot as hplot
from seispy.station import homestake, spiral
from seispy.station.stationdata import SeismometerArray
from seispy.station import StationArray
from scipy.linalg import pinv, svd
from matplotlib import rc
from seispy.test_mass import TestMass
import configparser
import optparse
import logging
rc('text', usetex=False)


def parse_command_line():
    """
    parse command line
    """
    parser = optparse.OptionParser()
    parser.add_option("--verbose", "-v",
                      help="verbosity", default=False,
                      dest="verbose", action="store_true")
    parser.add_option("--recovery-string", "-r",
                      help="recovery string", default='p',
                      dest="recovery_string", type=str)
    parser.add_option("--injection-file", "-I",
                      help="file with injection information", default=None,
                      dest="injection_file", type=str)
    parser.add_option("--inversion-condition", "-c",
                      help="smallest singular value to include", default=1e-4,
                      dest="inversion_condition", type=float)
    parser.add_option("--recovery-frequency", "-f",
                      help="recovery frequency", default=1,
                      dest="rec_freq", type=float)
    parser.add_option("--recovery-velocity-list", "-l",
                      help="comma separated list of velocities", default=5700,
                      dest="rec_velocity_list", type=str)
    parser.add_option("--nside",
                      help="nside for healpix", default=4,
                      dest="nside", type=str)
    parser.add_option("--make-recovery-plots",
                      help="include flag to make recovery plots", default=False,
                      dest="make_recovery_plots", action="store_true")
    parser.add_option("--duration",
                      help="duration of data to create", default=100,
                      dest="duration", type=float)
    parser.add_option("--save-tag",
                      help="prefix for plots", default='results',
                      dest="save_tag", type=str)

    params, args = parser.parse_args()
    return params


def run(params):
    # load up recovery velocities
    recovery_velocities = params.rec_velocity_list.strip().split(',')
    recovery_velocities = [float(r) for r in recovery_velocities]

    # logging
    logging.debug('Recovery Velocities:')
    logging.debug(recovery_velocities)
    logging.debug('Recovery string: %s', params.recovery_string)

    # error checking
    if len(recovery_velocities) != len(params.recovery_string):
        raise ValueError('N velocities should match N requested recoveries')

    # load injection/array file
    config = configparser.ConfigParser()
    print(params.injection_file)
    config.read(params.injection_file)

    # construct array for recovery
    if config['array']['array-type'] == 'spiral':
        stations = spiral(int(config['array']['nstations']),
                          radius=int(config['array']['radius']),
                          n_per_turn=int(config['array']['n_per_turn']))
    elif config['array']['array-type'] == 'homestake':
        try:
            station_names = config['array']['stations'].split(',')
            stations = homestake(station_names)
        except KeyError:
            stations = homestake()
    else:
        raise ValueError('Must specify an array-type "homestake" or "spiral"')

    # initialize data for array
    data = SeismometerArray.initialize_all_good(stations,
            float(config['array']['duration']),
                                                chans_type='fast_chans')

    # pop out array information

    # loop over what's left (injections)
    for key in config.keys():
        if key=='DEFAULT' or key=='array':
            continue
        if config[key]['type'] == 'p':
            data.add_p_wave(float(config[key]['amplitude']),
                            np.radians(float(config[key]['phi'])),
                            np.radians(float(config[key]['theta'])),
                            float(config['array']['frequency']),
                            float(config['array']['duration']),
                            float(config[key]['velocity']))
        if config[key]['type'] == 's':
            data.add_s_wave(config[key]['amplitude'],
                            np.radians(config[key]['phi']),
                            np.radians(config[key]['theta']),
                            config[key]['polarization'],
                            config[key]['frequency'],
                            config[key]['duration'],
                            c=config[key]['velocity'])
        if config[key]['type'] == 'r':
            data.add_r_wave(config[key]['amplitude'],
                            np.radians(config[key]['phi']),
                            np.radians(config[key]['theta']),
                            config[key]['frequency'],
                            config[key]['duration'],
                            c=config[key]['velocity'],
                            )
        if config[key]['type'] == 'l':
            data.add_l_wave(config[key]['amplitude'],
                            np.radians(config[key]['phi']),
                            np.radians(config[key]['theta']),
                            config[key]['frequency'],
                            config[key]['duration'],
                            c=config[key]['velocity'],
                            )
    # TODO: Auto add reflections if requested.
    csf = data.coherent_recovery(params.recovery_string,
                                 recovery_velocities,
                                 stations,
                                 frequency=float(config['array']['frequency']),
                                 nside=params.nside)
    freq = float(config['array']['frequency'])
    for rec in params.recovery_string:
        plt.figure(figsize=(12, 6))
        ax = plt.subplot(111, projection='mollweide')

        hplot.healpix_heatmap(np.abs(csf.maps[freq][rec]), cmap='plasma',
                              vmin=0, vmax=np.abs(csf.maps[freq][rec]).max())
        h = plt.colorbar()
        h.ax.set_ylabel('Amplitude of wave', fontsize=24)
        h.ax.tick_params(labelsize=16)
        plt.title("Amplitude Recovery", fontsize=32)
        hplot.outline_text(ax)
        ax.grid(True)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        ax.tick_params(labelsize=16)
        plt.savefig(params.save_tag+'-'+rec+'-amplitude.png')

    km_SI=1000
    csf.seismic['density'] = 2.5e3
    L = 40 * km_SI
    m = 40  # kg
    freqs = np.array([freq])
    x_ITMX = np.array([0, 0, 0])
    x_ITMY = np.array([0, 0, 0])
    x_ETMX = np.array([L, 0, 0])
    x_ETMY = np.array([0, L, 0])

    ITMX = TestMass('ITMX', x_ITMX, m, freqs)
    ITMY = TestMass('ITMY Underground', x_ITMY, m, freqs)
    ETMX = TestMass('ETMX', x_ETMX, m, freqs)
    ETMY = TestMass('ETMY', x_ETMY, m, freqs)

    # masses=[ITMX,ITMY,ETMX,ETMY]
    masses = [ITMX, ITMY]

    for mass in masses:
        print(mass.name)

        mass.get_acceleration_budget(csf)
        for ii, freq in enumerate(mass.freqs):
            a = mass.acceleration[freq]
            print('\t%2.2f Hz\n\t\tAx~ =%2.4e+i%2.4e m/s\n\t\tAy~ =%2.4e+i%2.4e m/s\n\t\tAz~ =%2.4e+i%2.4e m/s' % (freq, np.real(a[0]), np.imag(a[0]), np.real(a[1]), np.imag(a[1]), np.real(a[2]), np.imag(a[2])))
    newtonian_noise_sq = np.zeros(freqs.shape)
    newtonian_noise = np.zeros(freqs.shape)

    for ii, freq in enumerate(freqs):
        for mass in masses:
            a = mass.acceleration[freq]
            newtonian_noise_sq[ii] += np.abs(a[0] / L / (2 * np.pi * freq)**2)**2
        newtonian_noise[ii] = np.sqrt(newtonian_noise_sq[ii])

        print('%2.2f Hz\n\tNewtonian noise=%2.4e / sqrt(Hz)' %
              (freq, newtonian_noise[ii]))

    # reset to injection pixel and re-run NN
    for key in config.keys():
        if key == 'array' or key == 'DEFAULT':
            continue
        else:
            csf.maps[freq][config[key]['type']] = np.zeros(hp.nside2npix(params.nside))
    for key in config.keys():
        if key == 'array' or key == 'DEFAULT':
            continue
        else:
            phi, theta = np.radians(float(config[key]['phi'])), np.radians(float(config[key]['theta']))
            csf.maps[freq][config[key]['type']][hp.ang2pix(int(params.nside), theta, phi)] += float(config[key]['amplitude'])
    for rec in params.recovery_string:
        plt.figure(figsize=(12, 6))
        ax = plt.subplot(111, projection='mollweide')

        hplot.healpix_heatmap(np.abs(csf.maps[freq][rec]), cmap='plasma',
                              vmin=0, vmax=np.abs(csf.maps[freq][rec]).max())
        h = plt.colorbar()
        h.ax.set_ylabel('Amplitude of wave', fontsize=24)
        h.ax.tick_params(labelsize=16)
        plt.title("Amplitude Recovery", fontsize=32)
        hplot.outline_text(ax)
        ax.grid(True)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        ax.tick_params(labelsize=16)
        plt.savefig(params.save_tag+'-'+rec+'-injected-amplitude.png')


    for mass in masses:
        print(mass.name)

        mass.get_acceleration_budget(csf)
        for ii, freq in enumerate(mass.freqs):
            a = mass.acceleration[freq]
            print('\t%2.2f Hz\n\t\tAx~ =%2.4e+i%2.4e m/s\n\t\tAy~ =%2.4e+i%2.4e m/s\n\t\tAz~ =%2.4e+i%2.4e m/s' % (freq, np.real(a[0]), np.imag(a[0]), np.real(a[1]), np.imag(a[1]), np.real(a[2]), np.imag(a[2])))
    newtonian_noise_sq = np.zeros(freqs.shape)
    newtonian_noise = np.zeros(freqs.shape)

    for ii, freq in enumerate(freqs):
        for mass in masses:
            a = mass.acceleration[freq]
            newtonian_noise_sq[ii] += np.abs(a[0] / L / (2 * np.pi * freq)**2)**2
        newtonian_noise[ii] = np.sqrt(newtonian_noise_sq[ii])

        print('%2.2f Hz\n\tNewtonian noise=%2.4e / sqrt(Hz)' %
              (freq, newtonian_noise[ii]))



if __name__ == "__main__":
    params = parse_command_line()
    run(params)
