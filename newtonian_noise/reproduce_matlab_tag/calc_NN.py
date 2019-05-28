import numpy as np
from seispy.test_mass import TestMass
from seispy.station.coherent_field import CoherentSeismicField

G_SI=6.673e-11
km_SI=1000

freqs=np.array([0.1,1])
csf=CoherentSeismicField(freqs)

seismicfile = '../../matlab/anisotropic_NN/input/Rwave_params_20180423'
csf.readSeismicDataFromFile(seismicfile)

mapdir='../../HomestakeNewtonianNoiseEstimates/data_files/section_2/'
mapfile='coherent_map_p_s_reflected_waves_included.mat'
csf.readMapFromFile(mapdir,mapfile)

L=40*km_SI
m=40 #kg
x_ITMX=np.array([0,0,0])
x_ITMY=np.array([0,0,0])
x_ETMX=np.array([L,0,0])
x_ETMY=np.array([0,L,0])

ITMX=TestMass('ITMX',x_ITMX,m,freqs)
ITMY=TestMass('ITMY Underground',x_ITMY,m,freqs)
ETMX=TestMass('ETMX',x_ETMX,m,freqs)
ETMY=TestMass('ETMY',x_ETMY,m,freqs)

#masses=[ITMX,ITMY,ETMX,ETMY]
masses=[ITMX,ITMY]

for mass in masses:
    print(mass.name)

    mass.get_acceleration_budget(csf)
    for ii,freq in enumerate(mass.freqs):
        a=mass.acceleration[freq]
        print('\t%2.2f Hz\n\t\tAx~ =%2.4e+i%2.4e m/s\n\t\tAy~ =%2.4e+i%2.4e m/s\n\t\tAz~ =%2.4e+i%2.4e m/s'%(freq,
                                                                        np.real(a[0]),np.imag(a[0]),
                                                                         np.real(a[1]),np.imag(a[1]),
                                                                         np.real(a[2]),np.imag(a[2])))

newtonian_noise_sq=np.zeros(freqs.shape)
newtonian_noise=np.zeros(freqs.shape)

for ii,freq in enumerate(freqs):
    for mass in masses:
        a=mass.acceleration[freq]
        newtonian_noise_sq[ii]+=np.abs( a[0]/L/(2*np.pi*freq)**2 )**2
    newtonian_noise[ii]=np.sqrt(newtonian_noise_sq[ii])

    print('%2.2f Hz\n\tNewtonian noise=%2.4e / sqrt(Hz)'%(freq,newtonian_noise[ii]))

