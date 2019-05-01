import numpy as np
from scipy.io import loadmat

G_SI = 6.673e-11
km_SI = 1000


class CoherentSeismicField:
    """
    Data for seismic field
    Radiometer maps as a function of frequency and polarization
    Data about field: velocities, density, eigenfunctions, PSDs
    """

    def __init__(self, freqs):
        """
        Constructor
        """
        self.maps = {}  # Radiometer maps for different freqs and pols
        self.maps['thetas'] = None
        self.maps['phis'] = None

        self.asd = {}  # ASD (normalized to surface displacement)

        # Properties of seismic enviornment (velocities, eigenfunctions, PSDs)
        self.seismic = {}
        self.seismic['density'] = None

        self.polarizations = ['p', 'sh', 'sv', 'r']

        # Initialize frequency dependent properties of maps and seismic
        # if there is only one frequency past we need to make it
        # iterable in some way...do this here.
        self.freqs = freqs
        # check if numpy array
        if not isinstance(self.freqs, np.ndarray) and not isinstance(self.freqs, list):
            # if not, make it one
            self.freqs = np.array([self.freqs])
            # if shape is an empty tuple then it's still not
            # iterable so reshape it.
        for freq in self.freqs:
            self.maps[freq] = {}
            self.maps[freq]['p'] = None
            self.maps[freq]['sh'] = None
            self.maps[freq]['sv'] = None
            self.maps[freq]['r'] = None

            self.seismic[freq] = {}
            self.seismic[freq]['vp'] = None
            self.seismic[freq]['vs'] = None
            self.seismic[freq]['vr'] = None
            self.seismic[freq]['lambda'] = None
            self.seismic[freq]['mu'] = None
            # Should give H1,V1,H2,V2,h1,v1,h2,v2
            self.seismic[freq]['eigenfunction'] = None
            self.seismic[freq]['asd'] = None  # FUNCTION OF DEPTH??

            self.asd[freq] = None

    def setDefaultSeismicData(self):
        """
        Hardcoded data just to make life easier for testing...
        """
        pass

    def readSeismicDataFromFile(self, filename):
        """
        Read in R wave eignenfunction data
        Initialize other parameters using default values
        """

        feet = 0.3048  # meters

        # material properties
        self.seismic['density'] = 2.5e3  # % kg/m^3

        vs = 1500  # s wave velocity, m/s (TODO: update this)
        vp = 2300  # p wave velocity, m/s (TODO: Update this)

        # R-wave data
        c2, c4, Norm, a1, a2, a3, a4, lnVRfref, alpha = np.loadtxt(filename,
                                                                   skiprows=1, unpack=True)
        fref = 1

        for freq in self.freqs:
            self.seismic[freq]['vs'] = vs
            self.seismic[freq]['vp'] = vp
            self.seismic[freq]['vr'] = np.exp(lnVRfref) * (freq / fref)**alpha

            V1 = 1j * (1 - c2)
            V2 = 1j * c2
            H1 = Norm - c4
            H2 = c4
            krhomag = 2 * np.pi * freq / self.seismic[freq]['vr']
            v1 = krhomag * a1
            v2 = krhomag * a2
            h1 = krhomag * a3
            h2 = krhomag * a4
            self.seismic[freq]['eigenfunction'] = [
                H1, H2, V1, V2, h1, h2, v1, v2]

            #self.seismic[freq]['eigenfunction']['V1'] = 1j*(1-c2)
            #self.seismic[freq]['eigenfunction']['V2'] = 1j*c2
            #self.seismic[freq]['eigenfunction']['H1'] = Norm-c4
            #self.seismic[freq]['eigenfunction']['H2'] = c4

            #krhomag = 2*np.pi*freq/self.seismic[freq]['vr']
            #self.seismic[freq]['v1'] = krhomag*a1
            #self.seismic[freq]['v2'] = krhomag*a2
            #self.seismic[freq]['v3'] = krhomag*a3
            #self.seismic[freq]['v4'] = krhomag*a4

    def readMapFromFile_Incoherent_RemoveNegativePixels(self, mapdirs, mapfile):
        """
        Read map from file
        """

        for mapdir, freq in zip(mapdirs, self.freqs):
            # File name
            self.mapfile = mapdir + mapfile
            # Try to read the file
            try:
                tmp = loadmat(self.mapfile)
            except:
                raise Exception(
                    'CoherentSeismicField: NO SUCH FILE: %s' % mapfile)

            # Extract angles
            self.maps['thetas'] = tmp['thetas'][0]
            self.maps['phis'] = tmp['phis'][0]

            # Extract P-waves
            try:
                psd_map = tmp['maps']['p'][0][0][0]
                cut = psd_map < 0
                psd_map[cut] = 0
                self.maps[freq]['p'] = np.sqrt(psd_map)
                power = np.sum(np.abs(self.maps[freq]['p'])**2)
                #print('P-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO P waves')
                self.maps[freq]['p'] = None

            # Extract SH-waves
            try:
                psd_map = tmp['maps']['sh'][0][0][0]
                cut = psd_map < 0
                psd_map[cut] = 0
                self.maps[freq]['sh'] = np.sqrt(psd_map)
                power = np.sum(np.abs(self.maps[freq]['sh'])**2)
                #print('SH-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO SH waves')
                self.maps[freq]['sh'] = None

            # Extract SV-waves
            try:
                psd_map = tmp['maps']['sv'][0][0][0]
                cut = psd_map < 0
                psd_map[cut] = 0
                self.maps[freq]['sv'] = np.sqrt(psd_map)
                power = np.sum(np.abs(self.maps[freq]['sv'])**2)
                #print('SV-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO SV waves')
                self.maps[freq]['sv'] = None

            # Extract R-waves
            # ??? -- ONLY KEEP SURFACE?
            try:
                psd_map = tmp['maps']['r'][0][0][0]
                cut = psd_map < 0
                psd_map[cut] = 0
                self.maps[freq]['r'] = np.sqrt(psd_map)
                power = np.sum(np.abs(self.maps[freq]['r'])**2)

                cut = self.maps['thetas'] != np.pi / 2
                self.maps[freq]['r'][cut] = 0
                power2 = np.sum(np.abs(self.maps[freq]['r'])**2)
                norm = power / power2
                self.maps[freq]['r'] *= norm
                #print('R-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO R waves')
                self.maps[freq]['r'] = None

    def readMapFromFile(self, mapdir, mapfile):
        """
        Read map from file
        """

        # File name
        self.mapfile = mapdir + mapfile

        # Try to read the file
        try:
            tmp = loadmat(self.mapfile)
        except:
            raise Exception('CoherentSeismicField: NO SUCH FILE: %s' % mapfile)

        # Extract angles
        self.maps['thetas'] = tmp['thetas'][0]
        self.maps['phis'] = tmp['phis'][0]

        for freq in self.freqs:
            # Extract P-waves
            try:
                self.maps[freq]['p'] = tmp['maps']['p'][0][0][0]
                power = np.sum(np.abs(self.maps[freq]['p'])**2)
                #print('P-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO P waves')
                self.maps[freq]['p'] = None

            # Extract SH-waves
            try:
                self.maps[freq]['sh'] = tmp['maps']['sh'][0][0][0]
                power = np.sum(np.abs(self.maps[freq]['sh'])**2)
                #print('SH-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO SH waves')
                self.maps[freq]['sh'] = None

            # Extract SV-waves
            try:
                self.maps[freq]['sv'] = tmp['maps']['sv'][0][0][0]
                power = np.sum(np.abs(self.maps[freq]['sv'])**2)
                #print('SV-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO SV waves')
                self.maps[freq]['sv'] = None

            # Extract R-waves
            # ??? -- ONLY KEEP SURFACE?
            try:
                self.maps[freq]['r'] = tmp['maps']['r'][0][0][0]
                power = np.sum(np.abs(self.maps[freq]['r'])**2)

                cut = self.maps['thetas'] != np.pi / 2
                self.maps[freq]['r'][cut] = 0

                power2 = np.sum(np.abs(self.maps[freq]['r'])**2)

                norm = power / power2

                self.maps[freq]['r'][cut] *= norm

                #print('R-wave power: %e'%power)
            except:
                print('CoherentSeismicField Warning: NO R waves')
                self.maps[freq]['r'] = None

    def renormalizeMaps(self, asddirs, asdfilename):
        if len(asddirs) != len(self.freqs):
            raise Exception(
                'Length of asddirs should match the number of frequencies')
        for asddir, freq in zip(asddirs, self.freqs):
            self.read_total_station_power_from_file(asddir, asdfilename, freq)
        self.normalize_power()

    def read_total_station_power_from_file(self, asddir, asdfile, freq):
        if freq not in self.freqs:
            raise Exception('Given frequency not in CoherentSeismicField')
        tmp = loadmat(asddir + asdfile)
        total_station_power = tmp['total_station_power']
        self.asd[freq] = np.sqrt(np.max(total_station_power))

    def normalize_power(self):
        """
        Normalize the different polarizations, to a PSD at a given frequency
        Maintain relative power in each polarization
        For R waves, only keep surface
        """

        # loop over frequencies
        for freq in self.freqs:
            # sum power in each polarization for this frequency
            total_map_power_sq = 0
            for p in self.polarizations:
                total_map_power_sq += np.sum(np.abs(self.maps[freq][p])**2)

            # normalization constant needed so total power is the asd**2 at this frequency
            total_map_power = np.sqrt(total_map_power_sq)
            norm = self.asd[freq] / total_map_power

            # normalize each polarization at this frequency
            for p in self.polarizations:
                # sqrt(2) converts rms to p2p
                self.maps[freq][p] *= norm * np.sqrt(2)

    def get_amplitude(self, freq, pol):
        """
        Return a map of the amplitude, for a given polarization and frequency
        """
        return np.abs(self.maps[pol])

    def get_phase(self, freq, pol):
        """
        Return a map of the phase, for a given polarization and frequency
        """
        return np.angle(self.maps[freq][pol])
