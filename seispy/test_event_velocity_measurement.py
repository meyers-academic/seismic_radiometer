import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from seispy import event
from seispy.trace import Trace
import numpy as np
from astropy.table import Table, vstack
from scipy.stats import linregress
import statsmodels.formula.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import pandas as pd

tab = event.EventTable.read_seis_db('db/gary.db','db/pat.window')
pat = event.Event(*tab[9])
stations =\
    ['A4850','B4850','C4850','D4850','A2000','B2000','C2000','D2000','E2000',
    'YATES','ROSS','ORO','300','1700','A4100','C4100','D4100','RRDG']
tabs = []
envelopes = {}
data = {}
freqs = [0.4]
for station in stations:
    newtab, envelope, data1 =\
            pat.analyze(station,freqs,framedir='/archive/frames/homestake/',
            return_envelopes=True)
    tabs.append(newtab)
    envelopes[station]=envelope
    data[station]=data1
newtab = vstack(tabs)
filtered_table_HHR = newtab[newtab['channel']=='HHR']
print filtered_table_HHR
# plot velocity envelopes


for freq in freqs:
    for station in stations:
        dist = filtered_table_HHR[filtered_table_HHR['station']==station]['distance']
        vels = dist / (envelopes[station][freq]['HHR'].times.value -
                pat.time)
        plt.figure()
        plt.plot(vels, envelopes[station][freq]['HHR'].value)
        plt.savefig('env_vel_%s_%s' % (station, str(int(freq*10))))
        plt.close()
First = 1
st = envelopes['A4850'][0.4]['HHR'].times.value[0]
for key in envelopes.keys():
    if First:
        plot = (envelopes[key][0.4]['HHR']*1e9).plot(figsize=(6,12), alpha=0.5)
        plot.add_timeseries(data[key][0.4]['HHR']*1e9, alpha=0.5)
        ax = plot.gca()
        ax.set_xticks([])
        ax.set_xlabel('')
        First = 0
        ax.tick_params(
                axis='x',
                which='both',
                bottom='off',
                top='off',
                labelbottom='off')
        ax.tick_params(
                axis='y',
                labelsize=10,
                )
        ax.text(st+48, 0, key, fontsize=20)
        continue
    plot.add_timeseries((envelopes[key][0.4]['HHR'])*1e9, newax=True, alpha=0.5)
    plot.add_timeseries(data[key][0.4]['HHR']*1e9, color='k', alpha=0.5)
    ax = plot.gca()
    ax.set_xticks([])
    ax.set_xlabel('')
    ax.tick_params(
            axis='x',
            which='both',
            bottom='off',
            top='off',
            labelbottom='off')
    ax.tick_params(
            axis='y',
            labelsize=10,
            )
    ax.text(st+48, 0, key, fontsize=20)


plot.savefig('envelopes_test_HHR_short')
plot.close()
max_dist = np.max(filtered_table_HHR['distance'])
min_dist = np.max(filtered_table_HHR['distance'])
corr_idxs = []
delta_xs = []
time_lags = []
for station1 in stations:
    for station2 in stations:
        A2 =\
        np.hstack((np.zeros(envelopes[station2][0.4]['HHR'].size),envelopes[station2][0.4]['HHR'].value,
            np.zeros(envelopes[station2][0.4]['HHR'].size)))
        corr = np.correlate(A2,envelopes[station1][0.4]['HHR'].value)
        corr = corr[1:-1]
        times = envelopes[station2][0.4]['HHR'].times.value
        times -= times[0]
        times = times[1:]
        times = np.hstack((-times[::-1],[0],times[:]))
        corr_idxs.append(np.argmax(corr))
        time_lags.append(times[np.argmax(corr)])
print 'Done with getting time lags'
angles = np.arange(0,180,6)
vel_estimates = []
intercepts = []
stds = []
for angle in angles:
    delta_xs = []
    for station1 in stations:
        for station2 in stations:
            x1_vec = envelopes[station1][0.4]['HHR'].get_location()
            x2_vec = envelopes[station2][0.4]['HHR'].get_location()
            delta_x_vec = x1_vec - x2_vec
            omega = np.asarray([np.sin(np.radians(angle))*np.cos(np.radians(200)),
                np.sin(np.radians(angle))*np.sin(np.radians(200)),
                np.cos(np.radians(angle))])
            delta_x = np.dot(delta_x_vec, omega)
            delta_xs.append(delta_x)
    df = pd.DataFrame({'T':time_lags,'X':delta_xs})
    res = sm.ols(formula='X ~ T', data=df).fit()
    vel_estimates.append(res.params['T'])
    intercepts.append(res.params['Intercept'])
    print res.bse['T']
    stds.append(res.bse['T'])
print len(stds)

plt.figure()
plt.plot(angles, stds)
plt.xlabel('angle')
plt.ylabel('standard error')
plt.tight_layout()
plt.savefig('angles_best_fit')
plt.close()

plt.figure()
plt.plot(angles, vel_estimates)
plt.xlabel('angle')
plt.ylabel('velocity [m/s]')
plt.tight_layout()
plt.savefig('velocities_vs_angle')
plt.close()
stds = np.asarray(stds)
print stds.size
vel_estimates = np.asarray(vel_estimates)
#stds[np.where(vel_estimates < 0)] = np.inf

print 'Best fit angle:'
maxang = angles[np.argmin(stds)]
print maxang


print 'STDS argmin'
print np.argmin(stds)

delta_xs = []
heights = []
for station1 in stations:
    for station2 in stations:
        x1_vec = envelopes[station1][0.4]['HHR'].get_location()
        x2_vec = envelopes[station2][0.4]['HHR'].get_location()
        delta_x_vec = x1_vec - x2_vec
        print station1,station2,delta_x_vec
        heights.append(delta_x_vec[-1])
        omega = np.asarray([np.sin(np.radians(angle))*np.cos(np.radians(200)),
                np.sin(np.radians(angle))*np.sin(np.radians(200)),
                np.cos(np.radians(angle))])

        delta_x = np.dot(delta_x_vec, omega)
        delta_xs.append(delta_x)

time_lags = np.asarray(time_lags)
delta_xs = np.asarray(delta_xs)

df = pd.DataFrame({'T':time_lags,'X':delta_xs})
res = sm.ols(formula='X ~ T', data=df).fit()
prstd, iv_l, iv_u = wls_prediction_std(res)
print res.summary()

print 'Final Velocity Measurement:'
print res.conf_int()

plt.figure()
plt.scatter(time_lags, delta_xs, label='data', c=heights)
plt.plot(time_lags, res.predict(df['T']), label='best fit',
    color='r')
ax = plt.gca()
ax.plot(time_lags, iv_l, 'g--')
ax.plot(time_lags, iv_u, 'g--')
plt.legend(loc='best')
plt.ylabel('distance between stations [m]')
plt.xlabel('time delay [s]')
plt.title('Array prop. %4.2f degrees' % maxang)
plt.tight_layout()
plt.savefig('final_best_fit_distance_vs_time')
plt.close()

delta_xs = []
for station1 in stations:
    for station2 in stations:
        x1_vec = envelopes[station1][0.4]['HHR'].get_location()
        x2_vec = envelopes[station2][0.4]['HHR'].get_location()
        delta_x_vec = x1_vec - x2_vec
        omega = np.asarray([np.cos(np.radians(141)),
                np.sin(np.radians(141)),
                0])
        delta_x = np.dot(delta_x_vec, omega)
        delta_xs.append(delta_x)

time_lags = np.asarray(time_lags)
delta_xs = np.asarray(delta_xs)

df = pd.DataFrame({'T':time_lags,'X':delta_xs})
res = sm.ols(formula='X ~ T', data=df).fit()
prstd, iv_l, iv_u = wls_prediction_std(res)
print res.summary()

print '141 Velocity Measurement:'
print res.conf_int()

plt.figure()
plt.scatter(time_lags, delta_xs, label='data', c=heights)
plt.plot(time_lags, res.predict(df['T']), label='best fit',
    color='r')
ax = plt.gca()
ax.plot(time_lags, iv_l, 'g--')
ax.plot(time_lags, iv_u, 'g--')
plt.legend(loc='best')
plt.ylabel('distance between stations [m]')
plt.xlabel('time delay [s]')
plt.title('Array prop. %4.2f degrees' % 141)
plt.tight_layout()
plt.savefig('final_141_distance_vs_time')
plt.close()


delta_xs = []
for station1 in stations:
    for station2 in stations:
        x1_vec = envelopes[station1][0.4]['HHR'].get_location()
        x2_vec = envelopes[station2][0.4]['HHR'].get_location()
        delta_x_vec = x1_vec - x2_vec
        omega = np.asarray([np.sin(np.radians(135))*np.cos(np.radians(200)),
                np.sin(np.radians(135))*np.sin(np.radians(200)),
                np.cos(np.radians(135))])
        delta_x = np.dot(delta_x_vec, omega)
        delta_xs.append(delta_x)

time_lags = np.asarray(time_lags)
delta_xs = np.asarray(delta_xs)

df = pd.DataFrame({'T':time_lags,'X':delta_xs})
res = sm.ols(formula='X ~ T', data=df).fit()
prstd, iv_l, iv_u = wls_prediction_std(res)
print res.summary()

print '200 Velocity Measurement:'
print res.conf_int()

plt.figure()
plt.scatter(time_lags, delta_xs, label='data', c=heights)
plt.colorbar(label='vertical distance')
plt.plot(time_lags, res.predict(df['T']), label='best fit',
    color='r')
ax = plt.gca()
ax.plot(time_lags, iv_l, 'g--')
ax.plot(time_lags, iv_u, 'g--')
plt.legend(loc='best')
plt.ylabel('distance between stations [m]')
plt.xlabel('time delay [s]')
plt.title('Array prop. (200,135) degrees', fontsize=16)
plt.tight_layout()
plt.savefig('final_200_distance_vs_time')
plt.close()





# r-wave velocity
#First = 1
#st = envelopes[0.1]['HHR'].times.value[0]
#for key in envelopes.keys():
#    if First:
#        plot = (envelopes[key]['HHR']*1e9).plot(figsize=(6,12), alpha=0.5)
#        plot.add_timeseries(data[key]['HHR']*1e9, alpha=0.5)
#        ax = plot.gca()
#        ax.set_xticks([])
#        ax.set_xlabel('')
#        First = 0
#        ax.tick_params(
#                axis='x',
#                which='both',
#                bottom='off',
#                top='off',
#                labelbottom='off')
#        ax.tick_params(
#                axis='y',
#                labelsize=10,
#                )
#        ax.text(st+48, 0, key, fontsize=20)
#        continue
#    plot.add_timeseries((envelopes[key]['HHR'])*1e9, newax=True, alpha=0.5)
#    plot.add_timeseries(data[key]['HHR']*1e9, color='k', alpha=0.5)
#    ax = plot.gca()
#    ax.set_xticks([])
#    ax.set_xlabel('')
#    ax.tick_params(
#            axis='x',
#            which='both',
#            bottom='off',
#            top='off',
#            labelbottom='off')
#    ax.tick_params(
#            axis='y',
#            labelsize=10,
#            )
#    ax.text(st+48, 0, key, fontsize=20)
#
#plot.savefig('envelopes_test_HHR')
#plot.close()
