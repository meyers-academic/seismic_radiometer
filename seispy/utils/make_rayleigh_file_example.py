# script to show an example of making a rayleigh eigenfunction paramfile
# uses f=1Hz in table 6.7 of Tanner's thesis
import numpy as np

data={'C1':1,
      'C2':-0.05,
      'C3':1.95,
      'C4':-0.95,
      'a1':0.08,
      'a2':0.06,
      'a3':0.19,
      'a4':0.99,
      'v':1208}

np.save('rayleigh_paramfiles/default_rayleigh_params.npy',[data])
