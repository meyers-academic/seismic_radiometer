from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from seispy.station import SeismometerArray
from engine import parse_command_line
from seispy.station import homestake
import numpy as np
from collections import OrderedDict
import astropy.units as u
from seispy.seispy_io import print_params

st = 1116620022.63 +40
et = st + 100
rec_str = 'r'
v = [3000]
freq = 0.4
segdur = 25
nproc = 8



print 'Grabbing data...'
data = SeismometerArray.fetch_data(st, et,
        framedir='/archive/frames/homestake/', chans_type='fast_chans')

phimesh = float(15)
thetamesh = float(15)
thetas = np.arange(thetamesh, 180+thetamesh, thetamesh) * np.pi / 180
phis = np.arange(phimesh, 360+phimesh, phimesh) * np.pi / 180

stations = homestake()

# do recovery
print 'Doing recovery...'
maps, phis, thetas =\
        data.recovery_matrices(
            rec_str,
            stations,
            freq,
            v,
            fftlength=segdur,
            overlap=segdur/2.,
            autocorrelations=True,
            thetas=thetas, phis=phis,
            iter_lim=2000,
            nproc=nproc)
recovered_parameters = {}
for rec in maps.keys():
    recovered_parameters[rec] = OrderedDict()
    if rec == 'r':
        idx = np.where(thetas == np.pi / 2)
        final_power = maps[rec].data
        total_power =\
            (maps[rec].data/float(segdur)).sum()
        c = maps[rec].power_in_conf(0.5)
        c,p,t = maps[rec].get_contour(0.5)
        c = np.asarray(c).squeeze()
        recovered_parameters[rec]['phi: min, recovered, max']=np.array(p)*(180/np.pi)
        recovered_parameters[rec]['Total map power']=total_power * u.m**2
        recovered_parameters[rec]['recovered amplitude'] =\
            np.sqrt(np.max(final_power)/float(segdur)*u.m)
        plt.figure()
        plt.step(phis * (180 / np.pi), final_power, where='mid',label='mid')
        plt.scatter(phis * (180/np.pi), c*np.max(final_power))
        ax = plt.gca()
        ax.set_ylabel(r'Power [$\textrm{m}^2$]')
        ax.set_xlabel(r'$\phi$')
        plt.tight_layout()
        plt.savefig('%s/r_recovery_%s' %
                ('./','test_real_data'))
    else:
        maps[rec].data = maps[rec].data / float(segdur)
        tot_power = maps[rec].data.sum()
        p_rec = maps[rec].power_in_conf(0.2)
        contour, phi_vals, theta_vals = maps[rec].get_contour(0.20)
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
        plot.savefig('%s/%s_recovery_%s' % ('./', rec, 'test_real_data'))
        plot.close()
print_params(recovered_parameters)
