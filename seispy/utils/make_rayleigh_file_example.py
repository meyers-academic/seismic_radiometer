# script to show an example of making a rayleigh eigenfunction paramfile
# we use a toy example with arbitrarily chosen coefficients
import numpy as np

data={'C1':1,
      'C2':0.5,
      'C3':0.25,
      'C4':0.75,
      'a1':0.1,
      'a2':0.05,
      'a3':0.075,
      'a4':0.025,
      'v':1}

np.save('example_rayleigh_file.npy',[data])
