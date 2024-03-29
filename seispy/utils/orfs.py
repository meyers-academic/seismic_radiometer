from __future__ import division
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import sph_harm
from .utils import calc_travel_time2


def orf_p(ch1_vec, ch2_vec, det1_loc, det2_loc, vp, ff=None, thetamesh=1,
          phimesh=1):
    """
    Calculate p-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    ch2_vec : `list-like`
        channel2 vector
    det1_loc : `list-like`
        location of first sensor
    det2_loc : `list-like`
        location of second sensor
    vp : `float`
        velocity of p-wave

    Returns
    -------
    gamma : `numpy.ndarray`
        overlap reduction function for p-waves
    ff : `numpy.ndarray`
        frequency array
    """
    # get separation vector
    x_vec = np.array(det1_loc) - np.array(det2_loc)
    # make it a unit vector
    thetas = np.arange(0, 181, thetamesh) * np.pi / 180
    phis = np.arange(0, 360, phimesh) * np.pi / 180
    THETAS, PHIS = np.meshgrid(thetas, phis)
    dtheta = thetas[1] - thetas[0]
    dphi = phis[1] - phis[0]
    # much faster to vectorize things
    OmgX = np.sin(THETAS) * np.cos(PHIS)
    OmgY = np.sin(THETAS) * np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    if ff is None:
        ff = np.logspace(-2, 1, num=100)
    gammas = np.zeros(ff.size)
    for ii, f in enumerate(ff):
        gammas[ii] = np.sum(np.sum(np.sin(THETAS) *
                                   (OmgX * ch1_vec[0] +
                                    OmgY * ch1_vec[1] +
                                    OmgZ * ch1_vec[2]) *
                                   (OmgX * ch2_vec[0] +
                                    OmgY * ch2_vec[1] +
                                    OmgZ * ch2_vec[2]) *
                                   np.exp(2 * np.pi * 1j * f *
                                          (OmgX * x_vec[0] +
                                           OmgY * x_vec[1] +
                                           OmgZ * x_vec[0]) / vp) *
                                   (dtheta * dphi * 3 /
                                    (4 * np.pi))))
    return gammas, ff


def orf_s(ch1_vec, ch2_vec, det1_loc, det2_loc,
          vs, ff=None, thetamesh=1, phimesh=1):
    """
    Calculate p-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    ch2_vec : `list-like`
        channel2 vector
    det1_loc : `list-like`
        location of first sensor
    det2_loc : `list-like`
        location of second sensor
    vs : `float`
        velocity of s-wave

    Returns
    -------
    gamma : `numpy.ndarray`
        overlap reduction function for p-waves
    ff : `numpy.ndarray`
        frequency array
    """
    # get separation vector
    x_vec = det1_loc - det2_loc
    # make it a unit vector
    thetas = np.arange(0,181,thetamesh) * np.pi / 180
    phis = np.arange(0,360,phimesh) * np.pi / 180
    THETAS, PHIS = np.meshgrid(thetas, phis)
    dtheta = thetas[1] - thetas[0]
    dphi = phis[1] - phis[0]
    # much faster to vectorize things
    OmgX = np.sin(THETAS)*np.cos(PHIS)
    OmgY = np.sin(THETAS)*np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    # get vector perp to Omega
    PsiX = OmgX + 2 * np.random.rand(1)-1
    PsiY = OmgY + 2 * np.random.rand(1)-1
    # force dot product to be zero, make psi a unit vector
    PsiZ = -(PsiX*OmgX + PsiY*OmgY) / OmgZ
    normfac = np.sqrt(PsiX**2 + PsiY**2 + PsiZ**2)
    PsiX *= 1./normfac
    PsiY *= 1./normfac
    PsiZ *= 1./normfac
    # take cross product of omega and psi to
    # get third vector
    ChiX = OmgY*PsiZ - OmgZ*PsiY
    ChiY = OmgZ*PsiX - OmgX*PsiZ
    ChiZ = OmgX*PsiZ - OmgZ*PsiX

    if ff is None:
        ff = np.logspace(-2,1,num=100)
    gammas = np.zeros(ff.size)
    for ii, f in enumerate(ff):
        gammas[ii] = np.sum(np.sum(np.sin(THETAS) *
            (((PsiX * ch1_vec[0] +PsiY*ch1_vec[1] + PsiZ*ch1_vec[2]) *
            (PsiX * ch2_vec[0] +PsiY*ch2_vec[1] + PsiZ*ch2_vec[2])) +
            ((ChiX * ch1_vec[0] +ChiY*ch1_vec[1] + ChiZ*ch1_vec[2]) *
            (ChiX * ch2_vec[0] +ChiY*ch2_vec[1] + ChiZ*ch2_vec[2]))) *
            np.exp(2*np.pi*1j*f*(OmgX*x_vec[0] + OmgY*x_vec[1] +
            OmgZ*x_vec[0])/vs) * (dtheta*dphi*3/(8*np.pi))))
    return gammas, ff

def orf_p_directional2(ch1_vec, ch2_vec, det1_loc, det2_loc, vp, f, thetas=None,
                       phis=None, healpy=False):
    """
    Calculate p-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    ch2_vec : `list-like`
        channel2 vector
    det1_loc : `list-like`
        location of first sensor
    det2_loc : `list-like`
        location of second sensor
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
    Omg_shape = OmgX.shape
    s = OmgX.size
    OMEGA = np.zeros((s, 3))
    OMEGA[:, 0] = OmgX.flatten()
    OMEGA[:, 1] = OmgY.flatten()
    OMEGA[:, 2] = OmgZ.flatten()
    dt = calc_travel_time2(x_vec, OMEGA, vp)
    if dt.dtype == 'O':
        dt = np.array([dt[ii][0] for ii in range(dt.size)])
    dt = dt.reshape(Omg_shape)
    sf = ((OmgX * ch1_vec[0] +
           OmgY * ch1_vec[1] + OmgZ * ch1_vec[2]) *
          (OmgX * ch2_vec[0] +
           OmgY * ch2_vec[1] + OmgZ * ch2_vec[2]))
    gammas = sf * np.exp(2 * np.pi * 1j * f * dt)
    return gammas, phis, thetas


def orf_p_directional(ch1_vec, ch2_vec, det1_loc, det2_loc, vp, f, thetas=None,
                      phis=None, healpy=False):
    """
    Calculate p-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    ch2_vec : `list-like`
        channel2 vector
    det1_loc : `list-like`
        location of first sensor
    det2_loc : `list-like`
        location of second sensor
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
    Omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T
    dt = calc_travel_time2(x_vec, OMEGA, vp)
    if dt.dtype == 'O':
        dt = np.array([dt[ii][0] for ii in range(dt.size)])
    dt = dt.reshape(Omg_shape)
    sf = ((OmgX * ch1_vec[0] +
           OmgY * ch1_vec[1] + OmgZ * ch1_vec[2]) *
          (OmgX * ch2_vec[0] +
           OmgY * ch2_vec[1] + OmgZ * ch2_vec[2]))
    gammas = sf * np.exp(2 * np.pi * 1j * f * dt)
    return gammas, phis, thetas


def orf_s_directional(ch1_vec, ch2_vec, det1_loc, det2_loc, vs, f,
                      thetas=None, phis=None, healpy=False):
    """
    Calculate p-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    ch2_vec : `list-like`
        channel2 vector
    det1_loc : `list-like`
        location of first sensor
    det2_loc : `list-like`
        location of second sensor
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
    # dtheta = thetas[1] - thetas[0]
    # dphi = phis[1] - phis[0]
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
    sf1 = (PsiX * ch1_vec[0] + PsiY * ch1_vec[1] + PsiZ * ch1_vec[2]) *\
        (PsiX * ch2_vec[0] + PsiY * ch2_vec[1] + PsiZ * ch2_vec[2])
    sf2 = (ChiX * ch1_vec[0] + ChiY * ch1_vec[1] + ChiZ * ch1_vec[2]) *\
        (ChiX * ch2_vec[0] + ChiY * ch2_vec[1] + ChiZ * ch2_vec[2])
    omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T
    dt = calc_travel_time2(x_vec, OMEGA, vs).reshape(omg_shape)
    if dt.dtype == 'O':
        dt = np.array([dt[ii][0] for ii in range(dt.size)])
    gamma1 = sf1 * np.exp(2 * np.pi * 1j * f * dt)
    gamma2 = sf2 * np.exp(2 * np.pi * 1j * f * dt)
    return gamma1, gamma2, phis, thetas


def orf_l_directional(ch1_vec, ch2_vec, det1_loc, det2_loc, v, f,
                      thetas=None, phis=None, healpy=False,
                      decay_parameter=None):
    ch1_vec = ch1_vec / np.sqrt(np.dot(ch1_vec, ch1_vec))
    ch2_vec = ch2_vec / np.sqrt(np.dot(ch2_vec, ch2_vec))
    z1 = det1_loc[2]  # IS THIS CORRECT? (where is 'earth surface' z=0 defined)
    z2 = det2_loc[2]  # IS THIS CORRECT?
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
    sf1 = (PsiX * ch1_vec[0] + PsiY * ch1_vec[1] + PsiZ * ch1_vec[2]) * \
        (PsiX * ch2_vec[0] + PsiY * ch2_vec[1] + PsiZ * ch2_vec[2])
    # love waves should decay exponentially
    # and apply phase to gamma
    gamma = sf1 * np.exp(-decay_parameter * (k * z1 + k * z2)) * \
        np.exp(2 * np.pi * 1j * f * dt)
    return gamma, phis, thetas


def orf_r_directional(ch1_vec, ch2_vec, det1_loc, det2_loc, f,
                      thetas=None, phis=None, healpy=False,
                      rayleigh_paramfile=None,
                      rayleigh_paramdict=None, v=None):
    """
    Calculate r-wave overlap reduction function between
    two channels
    Parameters
    ----------
    ch1_vec : `list-like`
        channel 1 vector
    ch2_vec : `list-like`
        channel2 vector
    det1_loc : `list-like`
        location of first sensor
    det2_loc : `list-like`
        location of second sensor
    f : `float`
        frequency at which you would like the orf
    thetas: `list-like`
        OPTIONAL: list of theta coordinate values on which to evaluate ORF
    phis: `list-like`
        OPTIONAL: list of phi coordinates at which to evaluate ORF
    healpy: `boolean`
        OPTIONAL: use healpix? NOTE: if this is true,
        thetas and phis should also be given
    rayleigh_paramfile: `string`
        OPTIONAL: filename containing parameters for r wave eigenfunctions.
        If none is specified,
        default parameters will be used and a warning will be printed
    rayleigh_paramdict: `dict'
        OPTIONAL: dict containing values for r-wave eigenfunctions.
        parameters not specified will default to the values given
        in rayleigh_paramfile
    Returns
    -------
    gamma1 : `numpy.ndarray`
        overlap reduction function for s-waves, pol1
    phis : `numpy.ndarray`
        phis array
    thetas : `numpy.ndarray`
        theta values
    """
    # normalize channel 1 and 2 vectors
    ch1_vec = ch1_vec / np.sqrt(np.dot(ch1_vec, ch1_vec))
    ch2_vec = ch2_vec / np.sqrt(np.dot(ch2_vec, ch2_vec))

    # read parameters for eigenfunction
    # assumes bi-exponential model
    # given in equations 40,41 of seismic radiometer note 
    # https://zzz.physics.umn.edu/_media/groups/homestake/analysis/directional_analysis_v4.pdf
    # should be stored as a python dict, with values for these parameters
    # r1 = C1 exp(-a1*k*z) + C2 exp(-a2*k*z)
    # r2 = C3 exp(-a3*k*z) + C4 exp(-a4*k*z)
    # with k = omega/v = 2*pi*f/v
    if rayleigh_paramfile is None:
        #print 'WARNING: No Rayleigh paramfile specified, using default eigenfunction'
        a1 = 0.80
        a2 = 0.66
        a3 = 0.54
        a4 = 0.73
        if v is None:
            v = 2504
        C2 = -0.83
        C4 = -0.83
        Nvh = -0.78
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
        if v is None:
            v = data['v']
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
                v = rayleigh_paramdict['v']
    z1 = det1_loc[2]  # IS THIS CORRECT? (where is 'earth surface' z=0 defined)
    z2 = det2_loc[2]  # IS THIS CORRECT?
    k = 2 * np.pi * f / v

    r1_det1 = (np.exp(-a1 * k * z1) +
               C2 * np.exp(-a2 * k * z1)) * (1. / (1 + C2))
    r2_det1 = (np.exp(-a3 * k * z1) +
               C4 * np.exp(-a4 * k * z1)) * (Nvh / (1. + C4))
    r1_det2 = (np.exp(-a1 * k * z2) +
               C2 * np.exp(-a2 * k * z2)) * (1. / (1 + C2))
    r2_det2 = (np.exp(-a3 * k * z2) +
               C4 * np.exp(-a4 * k * z2)) * (Nvh / (1. + C4))
    # r1_det2=C1*np.exp(-a1*k*z2)+C2*np.exp(-a2*k*z2)
    # r2_det2=C3*np.exp(-a3*k*z2)+C4*np.exp(-a4*k*z2)
    # get separation vector
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

    omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T

    dt = calc_travel_time2(x_vec, OMEGA, v).reshape(omg_shape)
    if dt.dtype == 'O':
        dt = np.array([dt[ii][0] for ii in range(dt.size)])
    sf1 = (r1_det1 * np.dot(OMEGA, ch1_vec) +
           r2_det1 * np.exp(1j * np.pi / 2) * ch1_vec[2]).reshape(omg_shape)
    sf2 = (r1_det2 * np.dot(OMEGA, ch2_vec) +
           r2_det2 * np.exp(1j * np.pi / 2) * ch2_vec[2]).reshape(omg_shape)

    gamma = np.conj(sf1) * sf2 * np.exp(2 * np.pi * 1j * f * dt)

    return gamma, phis, thetas


def orf_picker(string, ch1_vec, ch2_vec,
               det1_loc, det2_loc, v, f, thetas=None,
               phis=None, healpy=False,
               rayleigh_paramfile=None,
               rayleigh_paramdict=None,
               decay_parameter=None
               ):
    if string == 'r':
        g1, p, t = orf_r_directional(ch1_vec, ch2_vec, det1_loc, det2_loc,
                                     f, thetas=thetas,
                                     phis=phis, healpy=healpy,
                                     rayleigh_paramfile=rayleigh_paramfile,
                                     rayleigh_paramdict=rayleigh_paramdict,
                                     v=v,
                                     )
        return g1.reshape((g1.size, 1)), g1.shape
    if string == 's':
        g1, g2, p, t = orf_s_directional(ch1_vec, ch2_vec, det1_loc, det2_loc,
                                         v, f, thetas=thetas, phis=phis,
                                         healpy=healpy)
        return g1.reshape((g1.size, 1)), g2.reshape((g2.size, 1)),\
            g1.shape, g2.shape
    if string == 'p':
        g1, p, t = orf_p_directional2(ch1_vec, ch2_vec, det1_loc, det2_loc,
                                      v, f, thetas=thetas,
                                      phis=phis, healpy=healpy)
        return g1.reshape((g1.size, 1)), g1.shape
    if string == 'l':
        g1, p, t = orf_l_directional(ch1_vec, ch2_vec, det1_loc, det2_loc,
                                     v, f, decay_parameter=decay_parameter, thetas=thetas,
                                     phis=phis, healpy=healpy)
        return g1.reshape((g1.size, 1)), g1.shape
