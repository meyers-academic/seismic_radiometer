import numpy as np 
from scipy.io import loadmat

G_SI=6.673e-11
km_SI=1000

class CoherentSeismicField:
    """
    Data for seismic field
    Radiometer maps as a function of frequency and polarization
    Data about field: velocities, density, eigenfunctions, PSDs
    """
    def __init__(self,freqs):
        """
        Constructor
        """
        self.maps={} # Radiometer maps for different frequencies and polarizations
        self.maps['thetas']=None
        self.maps['phis']=None

        self.seismic={} # Properties of seismic enviornment (velocities, eigenfunctions, PSDs)
        self.seismic['density']=None

        # Initialize frequency dependent properties of maps and seismic
        self.freqs=freqs
        for freq in self.freqs:
            self.maps[freq]={}
            self.maps[freq]['p']=None
            self.maps[freq]['sh']=None
            self.maps[freq]['sv']=None
            self.maps[freq]['r']=None

            self.seismic[freq]={}
            self.seismic[freq]['vp']=None
            self.seismic[freq]['vs']=None
            self.seismic[freq]['vr']=None
            self.seismic[freq]['lambda']=None
            self.seismic[freq]['mu']=None
            self.seismic[freq]['eigenfunction']=None # Should give H1,V1,H2,V2,h1,v1,h2,v2
            self.seismic[freq]['asd']=None # FUNCTION OF DEPTH??

    def setDefaultSeismicData(self):
        """
        Hardcoded data just to make life easier for testing...
        """
        pass

    def readSeismicDataFromFile(self,filename):
        """
        Read in R wave eignenfunction data
        Initialize other parameters using default values
        """

        feet=0.3048 # meters

        # material properties
        self.seismic['density']=2.5e3 # % kg/m^3


        vs=1500 # s wave velocity, m/s (TODO: update this)
        vp=2300 # p wave velocity, m/s (TODO: Update this)

        # R-wave data
        c2,c4,Norm,a1,a2,a3,a4,lnVRfref,alpha=np.loadtxt(filename,skiprows=1,unpack=True)
        fref=1

        for freq in self.freqs:
            self.seismic[freq]['vs']=vs
            self.seismic[freq]['vp']=vp
            self.seismic[freq]['vr']=np.exp(lnVRfref)*(freq/fref)**alpha

            self.seismic[freq]['V1']=1j*(1-c2)
            self.seismic[freq]['V2']=1j*c2
            self.seismic[freq]['H1']=Norm-c4
            self.seismic[freq]['H2']=c4

            krhomag=2*np.pi*freq/self.seismic[freq]['vr']
            self.seismic[freq]['v1']=krhomag*a1
            self.seismic[freq]['v2']=krhomag*a2
            self.seismic[freq]['v3']=krhomag*a3
            self.seismic[freq]['v4']=krhomag*a4


    def readMapFromFile(self,mapdir,mapfile):
        """
        Read map from file
        """

        # File name
        self.mapfile=mapdir+mapfile

        # Try to read the file
        try:
            tmp=loadmat(self.mapfile)
        except:
            raise Exception('CoherentSeismicField: NO SUCH FILE: %s'%mapfile)

        # Extract angles
        self.maps['thetas']=tmp['thetas'][0]
        self.maps['phis']=tmp['phis'][0]

        for freq in self.freqs:
            # Extract P-waves
            try:
                self.maps[freq]['p']=tmp['maps']['p'][0][0][0]
                power=np.sum(np.abs(self.maps[freq]['p'])**2)
                print('P-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO P waves')
                self.maps[freq]['p']=None

            # Extract SH-waves
            try:
                self.maps[freq]['sh']=tmp['maps']['sh'][0][0][0]
                power=np.sum(np.abs(self.maps[freq]['sh'])**2)
                print('SH-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO SH waves')
                self.maps[freq]['sh']=None

            # Extract SV-waves
            try:
                self.maps[freq]['sv']=tmp['maps']['sv'][0][0][0]
                power=np.sum(np.abs(self.maps[freq]['sv'])**2)
                print('SV-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO SV waves')
                self.maps[freq]['sv']=None

            # Extract R-waves
            # ??? -- ONLY KEEP SURFACE?
            try:
                self.maps[freq]['r']=tmp['maps']['r'][0][0][0]
                power=np.sum(np.abs(self.maps[freq]['r'])**2)
                print('R-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO R waves')
                self.maps[freq]['r']=None


    def normalize_power(self,depth):
        """
        Normalize the different polarizations, to a PSD at a given frequency and depth
        Maintain relative power in each polarization
        For R waves, only keep surface
        """
        pass

    def get_amplitude(self,freq,pol):
        """
        Return a map of the amplitude, for a given polarization and frequency
        """
        return np.abs(self.maps[pol])

    def get_phase(self,freq,pol):
        """
        Return a map of the phase, for a given polarization and frequency
        """
        return np.angle(self.maps[freq][pol])

