from __future__ import division
import numpy as np
from .utils import calc_travel_time2


def response_p_directional(ch1_vec, det1_loc, det2_loc, vp, f, thetas=None,
                           phis=None, healpy=False):
    """
    Calculate p-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    det1_loc : `list-like`
        location of sensor
    det2_loc : `list-like`
       origin location
    vp : `float`
        velocity of p-wave
    thetas : `numpy.ndarray`, optional
        list of theta values to use. Assumed to be angle FROM the north pole.
        If not supplied, defaults to 3 -> 177 in increments of 6 degrees.
    phis : `numpy.ndarray`, optional
        list of phi values to use. If not supplied,
        defaults to 3 -> 357 in increments of 6 degrees.

    Returns
    -------
    gamma : `numpy.ndarray`
        overlap reduction function for p-waves
    phis : `numpy.ndarray`
        phi values
    thetas : `numpy.ndarray`
        theta values
    """

    # get separation vector
    x_vec = + np.array(det1_loc) - np.array(det2_loc)
    # make it a unit vector
    if thetas is None:
        thetas = np.arange(3, 180, 6) * np.pi / 180
    if phis is None:
        phis = np.arange(3, 360, 6) * np.pi / 180
    if healpy:
        THETAS = thetas
        PHIS = phis
    else:
        THETAS, PHIS = np.meshgrid(thetas, phis)
    # much faster to vectorize things
    OmgX = np.sin(THETAS) * np.cos(PHIS)
    OmgY = np.sin(THETAS) * np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    Omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T
    dt = calc_travel_time2(x_vec, OMEGA, vp)
    if dt.dtype == 'O':
        dt = np.array([dt[ii][0] for ii in range(dt.size)])
    dt = dt.reshape(Omg_shape)
    sf = ((OmgX * ch1_vec[0] +
           OmgY * ch1_vec[1] + OmgZ * ch1_vec[2]))
    R = sf * np.exp(-2 * np.pi * 1j * f * dt)
    return R, phis, thetas


def response_s_directional(ch1_vec, det1_loc, det2_loc, vs, f,
                           thetas=None, phis=None, healpy=False):
    """
    Calculate p-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    det1_loc : `list-like`
        location of first sensor
    det2_loc : `list-like`
       origin location
    vs : `float`
        velocity of s-wave
    f : `float`
        frequency at which you would like the orf

    Returns
    -------
    gamma1 : `numpy.ndarray`
        overlap reduction function for s-waves, pol1
    gamma2 : `numpy.ndarray`
        overlap reduction function for s-waves, pol 2
    phis : `numpy.ndarray`
        phi values
    thetas : `numpy.ndarray`
        theta values
    """
    # get separation vector
    x_vec = np.array(det1_loc) - np.array(det2_loc)
    # make it a unit vector
    if thetas is None:
        thetas = np.arange(3, 180, 6) * np.pi / 180
    if phis is None:
        phis = np.arange(3, 360, 6) * np.pi / 180
    if healpy:
        THETAS = thetas
        PHIS = phis
    else:
        THETAS, PHIS = np.meshgrid(thetas, phis)
    # much faster to vectorize things
    OmgX = np.sin(THETAS) * np.cos(PHIS)
    OmgY = np.sin(THETAS) * np.sin(PHIS)
    OmgZ = (np.cos(THETAS))

    # get vector perp to Omega
    PsiX = -np.sin(PHIS)
    PsiY = np.cos(PHIS)

    # force dot product to be zero, make psi a unit vector
    PsiZ = np.zeros(PHIS.shape)

    # take cross product of omega and psi to
    # get third vector
    ChiX = -np.cos(THETAS) * np.cos(PHIS)
    ChiY = -np.cos(THETAS) * np.sin(PHIS)
    ChiZ = np.sin(THETAS)
    sf1 = (PsiX * ch1_vec[0] + PsiY * ch1_vec[1] + PsiZ * ch1_vec[2])
    sf2 = (ChiX * ch1_vec[0] + ChiY * ch1_vec[1] + ChiZ * ch1_vec[2])
    omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T
    dt = calc_travel_time2(x_vec, OMEGA, vs).reshape(omg_shape)
    if dt.dtype == 'O':
        dt = np.array([dt[ii][0] for ii in range(dt.size)])
    gamma1 = sf1 * np.exp(-2 * np.pi * 1j * f * dt)
    gamma2 = sf2 * np.exp(-2 * np.pi * 1j * f * dt)
    return gamma1, gamma2, phis, thetas


def response_l_directional(ch1_vec, det1_loc, det2_loc, v, f,
                           thetas=None, phis=None, healpy=False,
                           decay_parameter=0.26):
    """
    response matrix for l-wave.
    default decay parameter chosen for Homestake gold mine value
    estimated using both mine blasts and inversion from l-wave phase
    velocity estimates (i.e. to produce 1D s-wave velocity model)
    """
    ch1_vec = ch1_vec / np.sqrt(np.dot(ch1_vec, ch1_vec))
    z1 = det1_loc[2]  # IS THIS CORRECT? (where is 'earth surface' z=0 defined)
    k = 2 * np.pi * f / v

    x_vec = np.array(det1_loc) - np.array(det2_loc)
    # initialize angle arrays
    if thetas is None:
        thetas = np.arange(3, 180, 6) * np.pi / 180
    if phis is None:
        phis = np.arange(3, 360, 6) * np.pi / 180
    # define angle meshes. different for healpy vs non-healpy
    if healpy:
        THETAS = thetas
        PHIS = phis
    else:
        THETAS, PHIS = np.meshgrid(thetas, phis)
    # much faster to vectorize things
    # unit vector in direction of wave propagation
    # NOTE: OmgZ should be 0 for surface wave, but we'll allow it for now
    OmgX = np.sin(THETAS) * np.cos(PHIS)
    OmgY = np.sin(THETAS) * np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    # get vector perp to Omega and the surface
    # this is what we use for love waves (and Sh waves in another function)
    PsiX = -np.sin(PHIS)
    PsiY = np.cos(PHIS)
    # force dot product to be zero, make psi a unit vector
    PsiZ = np.zeros(PHIS.shape)

    omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T

    dt = calc_travel_time2(x_vec, OMEGA, v).reshape(omg_shape)
    if dt.dtype == 'O':
        dt = np.array([dt[ii][0] for ii in range(dt.size)])
    sf1 = (PsiX * ch1_vec[0] + PsiY * ch1_vec[1] + PsiZ * ch1_vec[2])
    # love waves should decay exponentially
    # and apply phase to gamma
    gamma = sf1 * np.exp(-decay_parameter * (k * z1)) * \
        np.exp(-2 * np.pi * 1j * f * dt)
    return gamma, phis, thetas


def response_r_directional(ch1_vec, det1_loc, det2_loc,
                           vr, f, rayleigh_paramdict=None,
                           rayleigh_paramfile=None,
                           thetas=None, phis=None, healpy=False):
    """
    Calculate r-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    det1_loc : `list-like`
        location of first sensor
    det2_loc : `list-like`
       origin location 
    vr : `float`
        velocity of s-wave
    f : `float`
        frequency at which you would like the orf

    Returns
    -------
    gamma1 : `numpy.ndarray`
        overlap reduction function for s-waves, pol1
    phis : `numpy.ndarray`
        phis array
    thetas : `numpy.ndarray`
        theta values
    """
    # get separation vector
    x_vec = np.array(det1_loc) - np.array(det2_loc)
    # make it a unit vector
    if thetas is None:
        thetas = np.arange(3, 180, 6) * np.pi / 180
    if phis is None:
        phis = np.arange(3, 360, 6) * np.pi / 180
    if healpy:
        THETAS = thetas
        PHIS = phis
    else:
        THETAS, PHIS = np.meshgrid(thetas, phis)
    # much faster to vectorize things
    OmgX = np.sin(THETAS) * np.cos(PHIS)
    OmgY = np.sin(THETAS) * np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    # only theta = pi/2 sticks around, okay? ok.
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T
    if rayleigh_paramfile is None:
        # print 'WARNING: No Rayleigh paramfile specified, using default eigenfunction'
        a1 = 0.47
        a2 = 0.73
        a3 = 1.51
        a4 = 0.25
        if vr is None:
            vr = 2504
        C2 = -1.29
        C4 = 2.29
        Nvh = -0.68
    else:
        data = np.load(rayleigh_paramfile)[0]
        # C1=data['C1']
        C2 = data['C2']
        # C3=data['C3']
        C4 = data['C4']
        a1 = data['a1']
        a2 = data['a2']
        a3 = data['a3']
        a4 = data['a4']
        if vr is None:
            vr = data['v']
    # if a parameter dict is specified, overwrite values
    if rayleigh_paramdict is None:
        pass
    else:
        for key in rayleigh_paramdict.keys():
            # if key=='C1':
            #     C1=rayleigh_paramdict['C1']
            if key == 'C2':
                C2 = rayleigh_paramdict['C2']
            # elif key=='C3':
            #     C3=rayleigh_paramdict['C3']
            elif key == 'C4':
                C4 = rayleigh_paramdict['C4']
            elif key == 'a1':
                a1 = rayleigh_paramdict['a1']
            elif key == 'a2':
                a2 = rayleigh_paramdict['a2']
            elif key == 'a3':
                a3 = rayleigh_paramdict['a3']
            elif key == 'a4':
                a4 = rayleigh_paramdict['a4']
            elif key == 'v':
                vr = rayleigh_paramdict['v']
    z1 = det1_loc[2]  # IS THIS CORRECT? (where is 'earth surface' z=0 defined)
    k = 2 * np.pi * f / vr

    r1_det1 = (np.exp(-a1 * k * z1) + C2 *
               np.exp(-a2 * k * z1)) * (1. / (1 + C2))
    r2_det1 = (np.exp(-a3 * k * z1) + C4 *
               np.exp(-a4 * k * z1)) * (Nvh / (1. + C4))

    # R-wave stuff
    # get vector perp to Omega
    # take cross product of omega and psi to
    # get third vector
    omg_shape = OmgX.shape
    sf1 = (r1_det1 * np.dot(OMEGA, ch1_vec) +
           r2_det1 * np.exp(1j * np.pi / 2) * ch1_vec[2]).reshape(omg_shape)
    dt = calc_travel_time2(x_vec, OMEGA, vr).reshape(omg_shape)
    if dt.dtype == 'O':
        dt = np.array([dt[ii][0] for ii in range(dt.size)])
    gamma = sf1 * np.exp(-2 * np.pi * 1j * f * dt)
    return gamma, phis, thetas


def response_picker(string, ch1_vec, det1_loc, origin_loc, velocity=None, frequency=None, thetas=None,
                    phis=None, epsilon=0.1, alpha=1000, healpy=False,
                    rayleigh_paramdict=None, rayleigh_paramfile=None,
                    decay_parameter=None):
    if velocity is None:
        raise ValueError('must specify velocity')
    if frequency is None:
        raise ValueError('must specify frequency')
    if string is 'r':
        g1, p, t = response_r_directional(ch1_vec, det1_loc, origin_loc,
                                          velocity, frequency,
                                          thetas=thetas, phis=phis,
                                          healpy=healpy,
                                          rayleigh_paramdict=rayleigh_paramdict,
                                          rayleigh_paramfile=rayleigh_paramfile)
        return g1.reshape((g1.size, 1)), g1.shape
    if string is 's':
        g1, g2, p, t = response_s_directional(ch1_vec, det1_loc, origin_loc,
                                              velocity, frequency,
                                              thetas=thetas,
                                              phis=phis, healpy=healpy)
        return (g1.reshape((g1.size, 1)), g2.reshape((g2.size, 1)), g1.shape,
                g2.shape)
    if string is 'p':
        g1, p, t = response_p_directional(ch1_vec, det1_loc, origin_loc,
                                          velocity, frequency,
                                          thetas=thetas, phis=phis,
                                          healpy=healpy)
        return g1.reshape((g1.size, 1)), g1.shape
    if string is 'l':
        g1, p, t = response_l_directional(ch1_vec, det1_loc, origin_loc,
                                          velocity, frequency,
                                          thetas=thetas,
                                          phis=phis, healpy=healpy,
                                          decay_parameter=decay_parameter)
        return g1.reshape((g1.size, 1)), g1.shape
