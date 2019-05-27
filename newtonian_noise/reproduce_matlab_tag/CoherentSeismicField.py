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

        self.asd={} # ASD (normalized to surface displacement)

        self.seismic={} # Properties of seismic enviornment (velocities, eigenfunctions, PSDs)
        self.seismic['density']=None

        self.polarizations=['p','sh','sv','r']

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

            self.asd[freq]=None


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

            V1=1j*(1-c2)
            V2=1j*c2
            H1=Norm-c4
            H2=c4
            krhomag=2*np.pi*freq/self.seismic[freq]['vr']
            v1=krhomag*a1
            v2=krhomag*a2
            h1=krhomag*a3
            h2=krhomag*a4
            self.seismic[freq]['eigenfunction']=[H1,H2,V1,V2,h1,h2,v1,v2]
            



    def readMapsFromFileList(self,mapdirs,mapfile,reg_method='coherent'):
        for mapdir,freq in zip(mapdirs,self.freqs):
            # File name
            mapfilename=mapdir+mapfile
            self.readMapFromFile(mapfilename,reg_method,freq)

    def readMapFromFile(self,mapfile,reg_method,freq):
        """
        Read map from file
        """

        # Try to read the file
        try:
            tmp=loadmat(mapfile)
        except:
            raise Exception('CoherentSeismicField: NO SUCH FILE: %s'%mapfile)


        # Extract angles
        self.maps['thetas']=tmp['thetas'][0]
        self.maps['phis']=tmp['phis'][0]

        for pol in self.polarizations:
            try:
                seismic_map=tmp['maps'][pol][0][0][0]
                regmap=self.regulateMap(seismic_map,pol,reg_method)
                self.maps[freq][pol]=regmap
            except:
                print('CoherentSeismicField Warning: NO %s waves'%pol)
                self.maps[freq][pol]=None

    def regulateMap(self,seismic_map,pol,method):
        """
        Given a seismic map, polarization, and a method, do any "cleaning" to raw map that is needed
        Note that incoherent maps are PSDs, and coherent maps give amplitudes, so we have to be careful about normalization

        INPUTS
        seismic_map: map from radiometer
        pol: polarization of map (surface waves are treated differently)
        method: if this is an incoherent map, how do we deal with negative pixels
        """
        def compute_power(seismic_map,is_coherent):
            if is_coherent:
                return np.sum(seismic_map)
            else:
                return np.sqrt(np.sum(seismic_map)**2)

        is_coherent= (method=='coherent')
        
        # if pol=='r':
        #     power1=compute_power(seismic_map,is_coherent)
        #     cut=self.maps['thetas']!=np.pi/2
        #     seismic_map[cut]=0
        #     power2=compute_power(seismic_map,is_coherent)
        #     norm=power1/power2
        #     seismic_map*=norm

        if method=='coherent':
            pass

        elif method=='Incoherent_LargestPixel':
            cut=seismic_map<max(seismic_map)
            seismic_map[cut]=0
            seismic_map=np.sqrt(seismic_map)


        elif method=='Incoherent_RemoveNegativePixels':
            cut=seismic_map<0
            seismic_map[cut]=0
            seismic_map=np.sqrt(seismic_map)

        elif method=='Incoherent_LargestPixel_Plus_Iso':
            midx=np.argmax(seismic_map)
            pmax=seismic_map[midx]
            #cut=seismic_map<0 # Matlab code does not remove negative pixels here
            #seismic_map[cut]=0
            pbcknd=compute_power(seismic_map,is_coherent)-pmax
            #if pol=='r':
            #    Npixels=np.sum(self.maps['thetas']!=np.pi/2)
            #else:
            #    Npixels=len(seismic_map)
            Npixels=len(seismic_map)
            seismic_map=np.ones(seismic_map.shape)*pbcknd/Npixels
            seismic_map[midx]+=pmax
            seismic_map=seismic_map+0*1j # lets complex square roots be taken appropriately
            seismic_map=np.sqrt(seismic_map) 
        return seismic_map


    def renormalizeMaps(self,asddirs,asdfilename):
        if len(asddirs) != len(self.freqs):
            raise Exception('Length of asddirs should match the number of frequencies')
        for asddir,freq in zip(asddirs,self.freqs):
            self.read_total_station_power_from_file(asddir,asdfilename,freq)
        

        self.normalize_surface_waves()
        self.normalize_power()



    def read_total_station_power_from_file(self,asddir,asdfile,freq):
        if freq not in self.freqs:
            raise Exception('Given frequency not in CoherentSeismicField')
        tmp=loadmat(asddir+asdfile)
        total_station_power=tmp['total_station_power']
        self.asd[freq]=np.sqrt(np.max(total_station_power))

    def normalize_surface_waves(self):
        # apply equator cut to r waves
        eq_cut=(self.maps['thetas']!=np.pi/2)
        for freq in self.freqs:
            rpower=np.sum(self.maps[freq]['r']**2)
            self.maps[freq]['r'][eq_cut]=0
            self.maps[freq]['r'] = self.maps[freq]['r'] * np.sqrt(rpower / np.sum(self.maps[freq]['r']**2))

    def normalize_power(self):
        """
        Normalize the different polarizations, to a PSD at a given frequency
        Maintain relative power in each polarization
        For R waves, only keep surface
        """

        # loop over frequencies
        for freq in self.freqs:
            # sum power in each polarization for this frequency
            total_map_power=self.asd[freq]**2
            map_power=0
            ppower={}
            for p in self.polarizations:
                ppower[p]=np.sum(self.maps[freq][p]**2) # no abs since we want to preserve minus signs to compare with matlab
                map_power+=ppower[p]
            norm=np.sqrt(total_map_power/map_power)
            # normalize each polarization at this frequency
            for p in self.polarizations:
                self.maps[freq][p] *= norm#* np.sqrt(2) # sqrt(2) converts rms to p2p

            print('F=%e Hz, seismic budget: Power=%e, R=%e, P=%e, SH=%e, SV=%e'%(freq,
                                                                                total_map_power,
                                                                                ppower['r']/total_map_power,
                                                                                ppower['p']/total_map_power,
                                                                                ppower['sh']/total_map_power,
                                                                                ppower['sv']/total_map_power))
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

