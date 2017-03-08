#! /usr/bin/python
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from seispy.station.stationdata import SeismometerArray
from seispy.station import StationArray, spiral, homestake
from seispy.seispy_io import read_config, print_params
import astropy.units as u
import os
import numpy as np
import ConfigParser
import optparse
from collections import OrderedDict

def parse_command_line():
    """
    parse command line
    """
    parser = optparse.OptionParser()
    parser.add_option("--config-file", "-c",
        help="configuration file", default=None,
        dest="config_file", type=str)
    params, args = parser.parse_args()
    return params

def main(params):
    print """
      _________      .__               .__         __________             .___.__                       __
     /   _____/ ____ |__| ______ _____ |__| ____   \______   \_____     __| _/|__| ____   _____   _____/  |_  ___________
     \_____  \_/ __ \|  |/  ___//     \|  |/ ___\   |       _/\__  \   / __ | |  |/  _ \ /     \_/ __ \   __\/ __ \_  __
     /        \  ___/|  |\___ \|  Y Y  \  \  \___   |    |   \ / __ \_/ /_/ | |  (  <_> )  Y Y  \  ___/|  | \  ___/|  | \/
    /_______  /\_____>__/______>__|_|__/__|\_____>  |____|___/(______/\_____| |__|\____/|__|_|__/\_____>__|  \_____>__|
    """
    print_params(params)
    try:
        os.mkdir(params['Recovery']['output_directory'])
    except:
        print 'Directory exists'
    if params['Recovery']['array_type']=='homestake':
        stations = homestake()
    else:
        stations = spiral(int(params['Recovery']['num_stations']))
    maxelevation=0
    # put things in terms of depth from highest station as opposed to
    # elevation above sea level
    for station in stations:
        if stations[station][2] > maxelevation:
            maxelevation = stations[station][2]
    for station in stations:
        stations[station][2] -= maxelevation
        stations[station][2] = (stations[station][2])
    # initialize
    data = SeismometerArray.initialize_all_good(stations,
        float(params['Recovery']['duration']))
    # do injections based on config file
    for key in params.keys():
        if key.split(' ')[0] == 'Injection':
            if params[key]['type']=='p':
                data.add_p_wave(float(params[key]['amplitude']),
                                float(params[key]['phi']) * np.pi / 180,
                                float(params[key]['theta']) * np.pi / 180,
                                float(params[key]['frequency']),
                                float(params['Recovery']['duration']),
                                c=float(params[key]['velocity']),
                                phase=float(params[key]['phase']),
                                Fs=float(params['Recovery']['sample_rate']))
            if params[key]['type']=='s':
                data.add_s_wave(float(params[key]['amplitude']),
                                float(params[key]['phi']) * np.pi / 180,
                                float(params[key]['theta']) * np.pi / 180,
                                float(params[key]['psi']) * np.pi / 180,
                                float(params[key]['frequency']),
                                float(params['Recovery']['duration']),
                                c=float(params[key]['velocity']),
                                phase=float(params[key]['phase']),
                                Fs=float(params['Recovery']['sample_rate']))
            if params[key]['type']=='r':
                data.add_r_wave(float(params[key]['amplitude']),
                                float(params[key]['phi']) * np.pi / 180,
                                float(params[key]['theta']) * np.pi / 180,
                                float(params[key]['epsilon']),
                                float(params[key]['alpha']),
                                float(params[key]['frequency']),
                                float(params['Recovery']['duration']),
                                c=float(params[key]['velocity']),
                                phase=float(params[key]['phase']),
                                Fs=float(params['Recovery']['sample_rate']))


    phimesh = float(params['Recovery']['phimesh'])
    thetamesh = float(params['Recovery']['phimesh'])
    thetas = np.arange(thetamesh, 180+thetamesh, thetamesh) * np.pi / 180
    phis = np.arange(phimesh, 360+phimesh, phimesh) * np.pi / 180

    # do recovery
    maps, phis, thetas =\
            data.recovery_matrices(
                params['Recovery']['recovery_string'],
                stations,
                params['Recovery']['frequency'],
                params['Recovery']['velocities'],
                fftlength=float(params['Recovery']['segdur']),
                overlap=float(params['Recovery']['segdur'])/2,
                autocorrelations=True,
                thetas=thetas, phis=phis,
                iter_lim=2000,
                nproc=int(params['Recovery']['nproc']),
                alpha=float(params['Recovery']['alpha']),
                epsilon=float(params['Recovery']['epsilon']))

    recovered_parameters = {}
    for rec in maps.keys():
        recovered_parameters[rec] = OrderedDict()
        if rec == 'r':
            idx = np.where(thetas == np.pi / 2)
            final_power = maps[rec].data
            total_power =\
                (maps[rec].data/float(params['Recovery']['segdur'])).sum()
            c = maps[rec].power_in_conf(0.5)
            c,p,t = maps[rec].get_contour(0.5)
            c = np.asarray(c).squeeze()
            recovered_parameters[rec]['phi: min, recovered, max']=np.array(p)*(180/np.pi)
            recovered_parameters[rec]['Total map power']=total_power * u.m**2
            recovered_parameters[rec]['recovered amplitude'] =\
                np.sqrt(np.max(final_power)/float(params['Recovery']['segdur'])*u.m)
            plt.figure()
            plt.step(phis * (180 / np.pi), final_power, where='mid',label='mid')
            plt.scatter(phis * (180/np.pi), c*np.max(final_power))
            ax = plt.gca()
            ax.set_ylabel(r'Power [$\textrm{m}^2$]')
            ax.set_xlabel(r'$\phi$')
            plt.savefig('%s/r_recovery_%s' %
                    (params['Recovery']['output_directory'],params['Recovery']['tag']))
        else:
            maps[rec].data = maps[rec].data / float(params['Recovery']['segdur'])
            tot_power = maps[rec].data.sum()
            p_rec = maps[rec].power_in_conf(0.5)
            contour, phi_vals, theta_vals = maps[rec].get_contour(0.50)
            recovered_parameters[rec]['phi low'] = phi_vals[0] * (180/np.pi)
            recovered_parameters[rec]['phi recovered'] = phi_vals[1] * (180/np.pi)
            recovered_parameters[rec]['phi max'] = phi_vals[2] * (180/np.pi)
            recovered_parameters[rec]['theta low'] = theta_vals[0]*(180/np.pi)
            recovered_parameters[rec]['theta recovered'] = theta_vals[1]*(180/np.pi)
            recovered_parameters[rec]['theta high'] = theta_vals[2]*(180/np.pi)
            recovered_parameters[rec]['Recovered amplitude'] = p_rec * u.m
            recovered_parameters[rec]['Total Map Power'] = tot_power * u.m**2
            plot = maps[rec].plot()
            ax = plot.gca()
            ax.set_title(r'$%s$-wave recovery' % rec, fontsize=12)
            plot.savefig('%s/%s_recovery_%s' % (params['Recovery']['output_directory'], rec, params['Recovery']['tag']))
            plot.close()
    print_params(recovered_parameters)

if __name__=="__main__":
    args = parse_command_line()
    params = read_config(args.config_file)
    main(params)

