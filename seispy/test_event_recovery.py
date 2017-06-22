import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from seispy import event
import numpy as np
from gwpy.table import Table


tab = event.EventTable.read_seis_db('db/gary.db','db/gary.window')
pat = event.Event(*tab[9])
newtab, envelopes, data =\
        pat.analyze('A4850',[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],framedir='/archive/frames/homestake/',
        return_envelopes=True)
filtered_HHR_table =  newtab['channel']=='HHR'
print filtered_HHR_table
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
#newtab.write('new_db',format='ascii')
#newtab2 = Table.read('new_db',format='ascii')
#print newtab2
