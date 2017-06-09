from __future__ import division
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import sph_harm
from .utils import calc_travel_time

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
        location of first sensor
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
    x_vec = np.array(det1_loc) - np.array(det2_loc)
    # make it a unit vector
    if thetas is None:
        thetas = np.arange(3,180,6) * np.pi / 180
    if phis is None:
        phis = np.arange(3,360,6) * np.pi / 180
    if healpy:
        THETAS = thetas
        PHIS = phis
    else:
        THETAS, PHIS = np.meshgrid(thetas, phis)
    # much faster to vectorize things
    OmgX = np.sin(THETAS)*np.cos(PHIS)
    OmgY = np.sin(THETAS)*np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    Omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T
    dt = calc_travel_time(x_vec, OMEGA, vp)
    dt = dt.reshape(Omg_shape)
    sf = ((OmgX*ch1_vec[0] +
            OmgY*ch1_vec[1] + OmgZ*ch1_vec[2]))
    R = sf *  np.exp(-2*np.pi*1j*f*dt)
    return R, phis, thetas

def response_s_directional(ch1_vec, det1_loc, det2_loc, vs, f,
        thetas=None,phis=None, healpy=False):
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
        thetas = np.arange(3,180,6) * np.pi / 180
    if phis is None:
        phis = np.arange(3,360,6) * np.pi / 180
    if healpy:
        THETAS = thetas
        PHIS = phis
    else:
        THETAS, PHIS = np.meshgrid(thetas, phis)
    dtheta = thetas[1] - thetas[0]
    dphi = phis[1] - phis[0]
    # much faster to vectorize things
    OmgX = np.sin(THETAS)*np.cos(PHIS)
    OmgY = np.sin(THETAS)*np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    # get vector perp to Omega
    PsiX = -np.sin(PHIS)
    PsiY = np.cos(PHIS)
    # force dot product to be zero, make psi a unit vector
    PsiZ = np.zeros(PHIS.shape)
    # take cross product of omega and psi to
    # get third vector
    ChiX = -np.cos(THETAS)*np.cos(PHIS)
    ChiY = -np.cos(THETAS)*np.sin(PHIS)
    ChiZ = np.sin(THETAS)
    sf1 = (PsiX*ch1_vec[0] +PsiY*ch1_vec[1] + PsiZ*ch1_vec[2])
    sf2 = (ChiX*ch1_vec[0] +ChiY*ch1_vec[1] + ChiZ*ch1_vec[2])
    omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T
    dt = calc_travel_time(x_vec, OMEGA, vs).reshape(omg_shape)
    gamma1 = sf1 * np.exp(-2*np.pi*1j*f*dt)
    gamma2 = sf2 * np.exp(-2*np.pi*1j*f*dt)
    return gamma1,gamma2,phis,thetas

def response_r_directional(ch1_vec, det1_loc, det2_loc, epsilon, alpha, vr, f,
        thetas=None,phis=None, healpy=False):
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
        thetas = np.arange(3,180,6) * np.pi / 180
    if phis is None:
        phis = np.arange(3,360,6) * np.pi / 180
    if healpy:
        THETAS = thetas
        PHIS = phis
    else:
        THETAS, PHIS = np.meshgrid(thetas, phis)
    # much faster to vectorize things
    OmgX = np.sin(THETAS)*np.cos(PHIS)
    OmgY = np.sin(THETAS)*np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    # only theta = pi/2 sticks around, okay? ok.
    OmgX[THETAS < np.pi / 2] = 0
    OmgY[THETAS < np.pi / 2] = 0
    OmgZ[THETAS < np.pi / 2] = 0
    OmgX[THETAS > np.pi / 2] = 0
    OmgY[THETAS > np.pi / 2] = 0
    OmgZ[THETAS > np.pi / 2] = 0
    # R-wave stuff
    R1 = np.cos(PHIS)
    R2 = np.sin(PHIS)
    R3 = np.exp(1j * np.pi / 2) * epsilon
    # get vector perp to Omega
    # take cross product of omega and psi to
    # get third vector
    sf1 = (R1*ch1_vec[0] + R2*ch1_vec[1] + R3*ch1_vec[2])

    omg_shape = OmgX.shape
    OMEGA = np.vstack((OmgX.flatten(), OmgY.flatten(), OmgZ.flatten())).T
    dt = calc_travel_time(x_vec, OMEGA, vr).reshape(omg_shape)
    gamma = sf1 * np.exp(-2*np.pi*1j*f*dt) * np.exp(-(det1_loc[2] +
        det2_loc[2]) / float(alpha))
    return gamma,phis,thetas

def response_picker(string, ch1_vec, det1_loc, origin_loc, v, f, thetas=None,
        phis=None, epsilon=0.1, alpha=1000, healpy=False):
    if string is 'r':
        g1, p, t =  response_r_directional(ch1_vec, det1_loc, origin_loc, epsilon,
                alpha, v, f, thetas=thetas, phis=phis, healpy=healpy)
        return g1.reshape((g1.size,1)), g1.shape
    if string is 's':
        g1, g2, p, t =  response_s_directional(ch1_vec, det1_loc, origin_loc,
                v, f, thetas=thetas, phis=phis, healpy=healpy)
        return g1.reshape((g1.size,1)), g2.reshape((g2.size,1)), g1.shape, g2.shape
    if string is 'p':
        g1, p, t = response_p_directional(ch1_vec, det1_loc, origin_loc,
                 v, f, thetas=thetas, phis=phis, healpy=healpy)
        return g1.reshape((g1.size,1)), g1.shape

