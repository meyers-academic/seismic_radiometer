import numpy as np
from TestMass import TestMass
from CoherentSeismicField import CoherentSeismicField
from matplotlib import pyplot as plt
from scipy.io import loadmat

G_SI=6.673e-11
km_SI=1000
feet = 0.3048 # meters

# 1117324816-1117332016
# 

matlab_results_dir='../../../matlab/anisotropic_NN/results/mats/'
filenames=['TIME_1117324816_DEPTH_0_REG_2_NN_budget.mat',
           'TIME_1117324816_DEPTH_800_REG_2_NN_budget.mat',
           'TIME_1117324816_DEPTH_4850_REG_2_NN_budget.mat']
matlab_results_files=[matlab_results_dir+x for x in filenames]
pol_idx_dict={'r':0,'p':1,'sh':2,'sv':3}
pol_color_dict={'r':'red','p':'blue','sh':'green','sv':'cyan'}
depth_color_dict={0:'red',800*feet:'blue',4850*feet:'green'}



freqs=np.arange(0.5,5.1,0.5)
csf=CoherentSeismicField(freqs)

seismicfile='../data/Rwave_params_20180423'
csf.readSeismicDataFromFile(seismicfile)

map_top_dir='../data/radiometer_maps/1117324816/'
map_data_dirs=['0_5_Hz/','1_0_Hz/','1_5_Hz/','2_0_Hz/','2_5_Hz/',
               '3_0_Hz/','3_5_Hz/','4_0_Hz/','4_5_Hz/','5_0_Hz/']
mapdirs=[map_top_dir+x for x in map_data_dirs]
mapfile='MAPS-1117324816-1117332016.mat'
#csf.readMapFromFile_Incoherent_RemoveNegativePixels(mapdirs,mapfile)
csf.readMapsFromFileList(mapdirs,mapfile,reg_method='Incoherent_LargestPixel_Plus_Iso')

asd_top_dir='./data/asds/1117324816/'
asd_data_dirs=['0_5_Hz/','1_0_Hz/','1_5_Hz/','2_0_Hz/','2_5_Hz/',
               '3_0_Hz/','3_5_Hz/','4_0_Hz/','4_5_Hz/','5_0_Hz/']
asd_dirs=[asd_top_dir+x for x in asd_data_dirs]
asd_filename='MAPS-1117324816-1117332016.mat'
csf.renormalizeMaps(asd_dirs,asd_filename)

L=40*km_SI
m=40 #kg
depths = feet*np.array([0,800,4850])

newtonian_noise={}
freqDet,ALIGO,ETD,CE,CEWB,CEPESS=np.loadtxt('data/CE_ASD_P1600143-v18.dat',unpack=True)

masses=[TestMass('ITMX',np.array([0,0,-depth]),m,freqs,is_coherent=False) for depth in depths]
newtonian_noise={depth:0 for depth in depths}
N_mirrors_per_arm=2
for mass,depth,matfile in zip(masses,depths,matlab_results_files):
    mass.get_acceleration_budget(csf)
    for p in mass.pols:
        newtonian_noise_sq=np.zeros(freqs.shape)
        for ii,f in enumerate(freqs):
            a=mass.acceleration_budget[f][p]
            newtonian_noise_sq[ii]+=(np.abs( a[0]/L/(2*np.pi*f)**2 )**2 
                                    +np.abs( a[1]/L/(2*np.pi*f)**2 )**2 )* N_mirrors_per_arm
        color=pol_color_dict[p]
        newtonian_noise_pol=np.sqrt(newtonian_noise_sq)
        newtonian_noise[depth]=np.sqrt(newtonian_noise[depth]**2+newtonian_noise_pol**2)
        plt.loglog(freqs,newtonian_noise_pol,'x--',label=p,color=color)
        tmp=loadmat(matfile)
        f_mat=tmp['freqs'][0]
        h_NN_pp=tmp['h_NN_pp']
        pol_idx=pol_idx_dict[p]
        plt.loglog(f_mat,h_NN_pp[:,pol_idx],':',label='MATLAB %s'%p,color=color)
    plt.loglog(freqs,newtonian_noise[depth],'o-',label='Total (%d ft)'%int(depth/feet),color='black')

    tmp=loadmat(matfile)
    f_mat=tmp['freqs'][0]
    h_NN=tmp['h_NN'][0]
    plt.loglog(f_mat,h_NN,':',label='MATLAB Total',color='black')




    plt.title(-mass.position[2])
    plt.loglog(freqDet,ETD,label='ETD')
    plt.loglog(freqDet,CE,label='CE')
    plt.xlim([0.5,10])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Strain/sqrt Hz')
    plt.grid()
    plt.legend()
    plt.ylim(1e-26,1e-20)
    plt.savefig('plots/pol_comparison_1117324816_depth_%2.4f.png'%depth)
    plt.clf()

for depth,matfile in zip(depths,matlab_results_files):
    color=depth_color_dict[depth]
    plt.loglog(freqs,newtonian_noise[depth],'o-',label='%d ft'%int(depth/feet),color=color)

    tmp=loadmat(matfile)
    f_mat=tmp['freqs'][0]
    h_NN=tmp['h_NN'][0]
    plt.loglog(f_mat,h_NN,':',label='MATLAB Total',color=color)
freqDet,ALIGO,ETD,CE,CEWB,CEPESS=np.loadtxt('data/CE_ASD_P1600143-v18.dat',unpack=True)
plt.loglog(freqDet,ETD,label='ETD')
plt.loglog(freqDet,CE,label='CE')
plt.xlim([0.5,10])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Strain/sqrt Hz')
plt.title('1117324816')
plt.grid()
plt.legend()
plt.savefig('plots/depth_comparison_1117324816.png')
