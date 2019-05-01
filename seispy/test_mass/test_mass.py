import numpy as np

G_SI = 6.673e-11
km_SI = 1000


class TestMass:
    def __init__(self, name, position, mass, freqs):
        """
        Constructor
        position -- 3d position vector x,y,z
        (z<0==>underground). Assume plane
            surface at =0
        mass -- kg
        freqs -- array of freqs in Hz at which to study test mass motion
        """
        self.name = name
        self.position = position
        self.mass = mass
        self.is_surface = (self.position[2] == 0)
        self.freqs = freqs
        self.acceleration = {}
        self.pols = ['p', 'r', 'sh', 'sv']
        self.acceleration_budget = {}
        for freq in self.freqs:
            self.acceleration[freq] = np.zeros(3)
            self.acceleration_budget[freq] = {}
            for pol in self.pols:
                self.acceleration_budget[freq][pol] = np.array([0, 0, 0])
        self.is_surface = (self.position[0] == 0)

    def get_acceleration(self):
        """
        Sum what is in the acceleration budget
        """
        for freq in self.freqs:
            for pol in self.pols:
                self.acceleration[freq] = (self.acceleration[freq] +
                                           self.acceleration_budget[freq][pol])

    def get_acceleration_budget(self, seismic_field):
        """
        Compute the acceleration budget of the test mass given a seismic field
        """
        maps = seismic_field.maps
        thetas = maps['thetas']
        phis = maps['phis']
        for freq in self.freqs:
            for pol in self.pols:
                skymap = maps[freq][pol]
                if skymap is None:
                    continue
                for theta, phi, a in zip(thetas, phis, skymap):
                    self.acceleration_budget[freq][pol] = (self.acceleration_budget[freq][pol] +
                                                           self.get_acceleration_component(freq, theta, phi, pol, a, seismic_field))
                self.acceleration[freq] = (self.acceleration[freq] +
                                           self.acceleration_budget[freq][pol])

    def get_acceleration_component(self, freq, theta, phi, pol, a, seismic_field):
        """
        Compute the acceleration component
        """
        if pol == 'r':
            if self.is_surface:
                return self.get_acceleration_R_surface(freq, theta, phi, a, seismic_field.seismic)
            else:
                return self.get_acceleration_R_cavity(freq, theta, phi, a, seismic_field.seismic)
        elif pol == 'p':
            if self.is_surface:
                return self.get_acceleration_P_surface(freq, theta, phi, a, seismic_field.seismic)
            else:
                return self.get_acceleration_P_cavity(freq, theta, phi, a, seismic_field.seismic)
        elif pol == 'sh':
            if self.is_surface:
                return self.get_acceleration_SH_surface(freq, theta, phi, a, seismic_field.seismic)
            else:
                return self.get_acceleration_SH_cavity(freq, theta, phi, a, seismic_field.seismic)
        elif pol == 'sv':
            if self.is_surface:
                return self.get_acceleration_SV_surface(freq, theta, phi, a, seismic_field.seismic)
            else:
                return self.get_acceleration_SV_cavity(freq, theta, phi, a, seismic_field.seismic)

        raise Exception('ILLEGAL POLARIZATION: %s' % pol)

    def get_acceleration_R_surface(self, freq, theta, phi, a, seismic_field):
        """
        Flag 9 in GGN.m
        incident Rayleigh wave, using Jan's equation for above the surface

        FT convention: int dt e^i*omega*t h(t)
        Can see by comparing the t terms in GGN
        """
        rho = self.position * \
            np.array(
                [1, 1, 0])  # projection of the position on the horizontal plane
        h = -self.position[2]  # depth of the point

        dens = seismic_field['density']

        speed = seismic_field[freq]['vr']
        wavelength = speed / freq
        k = 2 * np.pi / wavelength * np.array([np.sin(theta) * np.sin(phi),
                                               np.sin(theta) * np.cos(phi),
                                               np.cos(theta)])

        krho = k * np.array([1, 1, 0])  # horizontal wavenumber
        # magnitude of the horizontal wavenumber
        krhomag = np.linalg.norm(krho)
        kmag = np.linalg.norm(k)  # magnitude of the total wavenumber

        # THIS INFORMATION SHOULD BE STORED IN
        # seismic[freq]['eigenfunction']
        # seismic[freq]['vr']
        H1, H2, V1, V2, h1, h2, v1, v2 = seismic_field[freq]['eigenfunction']
        AH1 = H1 * a
        AH2 = H2 * a
        AV1 = V1 * a
        AV2 = V2 * a
        #AH1 = A(1)
        #AH2 = A(2)
        #AV1 = A(3)
        #AV2 = A(4)
        #h1 = A(5)
        #h2 = A(6)
        #v1 = A(7)
        #v2 = A(8)
        #speed = A(9)
        omega = speed * kmag

        del1 = (-2 * np.pi * G_SI * dens *
                np.exp(-np.dot(krhomag, h)) *
                np.exp(1j * np.dot(krho, rho)))
        del2 = (AH1 / (h1 + krhomag) +
                AH2 / (h2 + krhomag) +
                AV1 / (v1 + krhomag) +
                AV2 / (v2 + krhomag))
        delphi = del1 * del2
        ax = -1j * krho[0] * delphi
        ay = -1j * krho[1] * delphi
        az = -krhomag * delphi
        dela = np.array([ax, ay, az])

        return dela

    def get_acceleration_R_cavity(self, freq, theta, phi, a, seismic_field):
        """
        Flags 7,8 in GGN.m
        incident Rayleigh wave, using Jan's equation for BELOW the surface
        incident Rayleigh wave interacting with the surface of a cavity
        """
        rho = self.position * \
            np.array(
                [1, 1, 0])  # projection of the position on the horizontal plane
        h = -self.position[2]  # depth of the point
        x = self.position

        dens = seismic_field['density']

        speed = seismic_field[freq]['vr']
        wavelength = speed / freq
        k = 2 * np.pi / wavelength * np.array([np.sin(theta) * np.sin(phi),
                                               np.sin(theta) * np.cos(phi),
                                               np.cos(theta)])

        krho = k * np.array([1, 1, 0])  # horizontal wavenumber
        # magnitude of the horizontal wavenumber
        krhomag = np.linalg.norm(krho)
        kmag = np.linalg.norm(k)  # magnitude of the total wavenumber

        # THIS INFORMATION SHOULD BE STORED IN
        # seismic[freq]['eigenfunction']
        # seismic[freq]['vr']
        H1, H2, V1, V2, h1, h2, v1, v2 = seismic_field[freq]['eigenfunction']
        AH1 = H1 * a
        AH2 = H2 * a
        AV1 = V1 * a
        AV2 = V2 * a
        omega = speed * kmag

        # incident Rayleigh wave, using Jan's equation for BELOW the surface
        del1 = 2 * np.pi * G_SI * dens * np.exp(1j * np.dot(krho, rho))
        del2 = np.exp(-h * krhomag) * (AH1 / (-h1 + krhomag) +
                                       AH2 / (-h2 + krhomag) -
                                       AV1 / (-v1 + krhomag) -
                                       AV2 / (-v2 + krhomag))
        del3 = (2 * AH1 * krhomag * np.exp(-h1 * h) / (h1**2 - krhomag**2) +
                2 * AH2 * krhomag * np.exp(-h2 * h) / (h2**2 - krhomag**2) +
                2 * AV1 * v1 * np.exp(-v1 * h) / (krhomag**2 - v1**2) +
                2 * AV2 * v2 * np.exp(-v2 * h) / (krhomag**2 - v2**2))
        del4 = (2 * h1 * AH1 * krhomag * np.exp(-h1 * h) / (h1**2 - krhomag**2) +
                2 * h2 * AH2 * krhomag * np.exp(-h2 * h) / (h2**2 - krhomag**2) +
                2 * v1 * AV1 * v1 * np.exp(-v1 * h) / (krhomag**2 - v1**2) +
                2 * v2 * AV2 * v2 * np.exp(-v2 * h) / (krhomag**2 - v2**2))
        ax = -1j * krho[0] * del1 * (del3 + del2)
        ay = -1j * krho[1] * del1 * (del3 + del2)
        az = -del1 * (del4 + krhomag * del2)
        dela = np.array([ax, ay, az])

        # incident Rayleigh wave interacting with the surface of a cavity
        del5 = -4 * np.pi * G_SI * dens / 3 * np.exp(1j * (np.dot(k, x)))
        ax = del5 * 1j * (AH1 * np.exp(-h1 * h) + AH2 *
                          np.exp(-h2 * h)) * k[0] / kmag
        ay = del5 * 1j * (AH1 * np.exp(-h1 * h) + AH2 *
                          np.exp(-h2 * h)) * k[1] / kmag
        az = del5 * (AV1 * np.exp(-v1 * h) + AV2 * np.exp(-v2 * h))
        dela = dela + np.array([ax, ay, az])

        return dela

    def get_acceleration_P_surface(self, freq, theta, phi, a, seismic_field):
        """
        Flags 11,13 in GGN.m
        p wave, infinite space, bulk contribution. eqn 52 of 1507.05850 (see also text after eqn 61)
        p wave, half space, reflection from surface. surface term of eqn 84 of 1507.05850 
        """
        x = self.position
        dens = seismic_field['density']
        h = -self.position[2]  # depth of the point

        speed = seismic_field[freq]['vp']
        wavelength = speed / freq
        k = 2 * np.pi / wavelength * np.array([np.sin(theta) * np.sin(phi),
                                               np.sin(theta) * np.cos(phi),
                                               np.cos(theta)])
        krho = k * np.array([1, 1, 0])  # horizontal wavenumber
        # magnitude of the horizontal wavenumber
        krhomag = np.linalg.norm(krho)

        kmag = np.linalg.norm(k)  # magnitude of the total wavenumber

        # bulk
        nhat = k / kmag  # dispacement vector
        ax = 4 * np.pi * G_SI * dens * a * nhat[0] * np.exp(1j * np.dot(k, x))
        ay = 4 * np.pi * G_SI * dens * a * nhat[1] * np.exp(1j * np.dot(k, x))
        az = 4 * np.pi * G_SI * dens * a * nhat[2] * np.exp(1j * np.dot(k, x))
        dela = np.array([ax, ay, az])

        # surface
        delphi = (2 * np.pi * G_SI * dens * a * 1 / (1j * kmag) * np.exp(-krhomag * h) *
                  np.exp(np.dot(1j * krho, self.position)))
        ax = -1j * krho[0] * delphi
        ay = -1j * krho[1] * delphi
        az = -krhomag * delphi  # d/dz = - d/dh
        dela = dela + np.array([ax, ay, az])

        return dela

    def get_acceleration_P_cavity(self, freq, theta, phi, a, seismic_field):
        """
        Flags 11,12,13 in GGN.m
        p wave, infinite space, bulk contribution. eqn 52 of 1507.05850 (see also text after eqn 61)
        p wave, infinite space, reflection from small cavity centered at x. eqn 61 of 1507.05850
        p wave, half space, reflection from surface. surface term of eqn 84 of 1507.05850 
        """

        # Flags 11,13 already implemented in the surface part
        dela = self.get_acceleration_P_surface(
            freq, theta, phi, a, seismic_field)

        x = self.position
        speed = seismic_field[freq]['vp']
        wavelength = speed / freq
        dens = seismic_field['density']
        k = 2 * np.pi / wavelength * np.array([np.sin(theta) * np.sin(phi),
                                               np.sin(theta) * np.cos(phi),
                                               np.cos(theta)])

        kmag = np.linalg.norm(k)  # magnitude of the total wavenumber
        nhat = k / kmag  # dispacement vector

        ax = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[0] * np.exp(1j * np.dot(k, x))
        ay = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[1] * np.exp(1j * np.dot(k, x))
        az = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[2] * np.exp(1j * np.dot(k, x))
        dela = dela + np.array([ax, ay, az])

        return dela

    def get_acceleration_SH_surface(self, freq, theta, phi, a, seismic_field):
        """
        Flag 17 in GGN.m
        surface contribution of sh wave
        """
        dela = np.zeros(3)
        return dela

    def get_acceleration_SH_cavity(self, freq, theta, phi, a, seismic_field):
        """
        Flag 14 in GGN.m
        """
        dens = seismic_field['density']

        x = self.position
        speed = seismic_field[freq]['vs']
        wavelength = speed / freq
        k = 2 * np.pi / wavelength * np.array([np.sin(theta) * np.sin(phi),
                                               np.sin(theta) * np.cos(phi),
                                               np.cos(theta)])

        nhat = np.array([-k[2], k[1], 0])  # dispacement vector
        nhat = nhat / np.linalg.norm(nhat)

        ax = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[0] * np.exp(1j * np.dot(k, x))
        ay = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[1] * np.exp(1j * np.dot(k, x))
        az = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[2] * np.exp(1j * np.dot(k, x))
        dela = np.array([ax, ay, az])

        return dela

    def get_acceleration_SV_surface(self, freq, theta, phi, a, seismic_field):
        """
        Flags 16 in GGN.m
        sv wave, half space, reflection from surface. eqn 85 of 1507.05850

        NOTE -- Flags 13, 16 are the same (except speeds), double check this is right
        """
        dens = seismic_field['density']
        speed = seismic_field[freq]['vs']
        wavelength = speed / freq
        k = 2 * np.pi / wavelength * np.array([np.sin(theta) * np.sin(phi),
                                               np.sin(theta) * np.cos(phi),
                                               np.cos(theta)])
        kmag = np.linalg.norm(k)
        krho = k * np.array([1, 1, 0])  # horizontal wavenumber
        # magnitude of the horizontal wavenumber
        krhomag = np.linalg.norm(krho)
        h = -self.position[2]  # depth of the point

        delphi = (2 * np.pi * G_SI * dens * a * 1 / (1j * kmag) * np.exp(-krhomag * h) *
                  np.exp(np.dot(1j * krho, self.position)))
        ax = -1j * krho[0] * delphi
        ay = -1j * krho[1] * delphi
        az = -krhomag * delphi  # d/dz = - d/dh
        dela = np.array([ax, ay, az])

        return dela

    def get_acceleration_SV_cavity(self, freq, theta, phi, a, seismic_field):
        """
        Flags 15,16 in GGN.m
        sv wave, half space, reflection from surface. eqn 85 of 1507.05850
        sv wave, infinite space, reflection from small cavity centered at x. eqn 62 of 1507.05850

        NOTE -- Flags 14, 15 are the same(except speeds), double check this is right
        """
        # Flag 16 already implemented in the surface part
        dela = self.get_acceleration_SV_surface(
            freq, theta, phi, a, seismic_field)

        x = self.position
        dens = seismic_field['density']
        h = -self.position[2]  # depth of the point

        speed = seismic_field[freq]['vs']
        wavelength = speed / freq
        k = 2 * np.pi / wavelength * np.array([np.sin(theta) * np.sin(phi),
                                               np.sin(theta) * np.cos(phi),
                                               np.cos(theta)])
        # krho = k * np.array([1, 1, 0])  #horizontal wavenumber
        # krhomag = np.linalg.norm(krho)  #magnitude of the horizontal wavenumber

        # kmag = np.linalg.norm(k)  #magnitude of the total wavenumber

        nhat = np.array([-k[2], k[1], 0])  # dispacement vector
        nhat = nhat / np.linalg.norm(nhat)

        ax = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[0] * np.exp(1j * np.dot(k, x))
        ay = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[1] * np.exp(1j * np.dot(k, x))
        az = -4 / 3 * np.pi * G_SI * dens * a * \
            nhat[2] * np.exp(1j * np.dot(k, x))
        dela = dela + np.array([ax, ay, az])

        return dela

    def shift_position_parallel_to_surface(self, shiftx, shifty, seismic_field):
        """
        Shift the position relative to surface
        """
        self.position[0] += shiftx
        self.position[1] += shifty
        # Need to loop over each component in acceleration (freq,pol)
        # Compute the wavenumber for that polarization,freq
        # Multiply acceleration component by a phase shift e^(i k.x)
        # Resum components to get total accleration
        # OBSTRUCTION... k.x depends on sky position...
        # do we need to store the budget for each pixel too?
        v = {}
        v['sh'] = seismic_field['vs']
        v['sh'] = seismic_field['vs']
        v['p'] = seismic_field['vp']
        v['r'] = seismic_field['vr']
        for freq in self.freqs:
            for pol in self.pols:
                k = 2 * np.pi * v[pol] / freq
                self.acceleration_budget[freq][pol] *= 1
        self.get_acceleration()

    def print_acceleration_to_file(self, filename):
        """
        Do we actually need this?
        """
        pass
