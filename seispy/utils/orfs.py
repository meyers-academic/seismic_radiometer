from __future__ import division
import numpy as np
import time
from scipy.interpolate import interp1d
from scipy.special import sph_harm

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
    if ff is None:
        ff = np.logspace(-2,1,num=100)
    gammas = np.zeros(ff.size)
    for ii,f in enumerate(ff):
        gammas[ii] = np.sum(np.sum(np.sin(THETAS) * (OmgX*ch1_vec[0] +
            OmgY*ch1_vec[1] + OmgZ*ch1_vec[2]) * (OmgX*ch2_vec[0] +
            OmgY*ch2_vec[1] + OmgZ*ch2_vec[2]) *
            np.exp(2*np.pi*1j*f*(OmgX*x_vec[0] + OmgY*x_vec[1] +
                OmgZ*x_vec[0])/vp) * (dtheta*dphi*3/(4*np.pi))))
    return gammas, ff

def orf_s(ch1_vec, ch2_vec, det1_loc, det2_loc, vs, ff=None, thetamesh=1,phimesh=1):
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
    for ii,f in enumerate(ff):
        gammas[ii] = np.sum(np.sum(np.sin(THETAS) *
            (((PsiX*ch1_vec[0] +PsiY*ch1_vec[1] + PsiZ*ch1_vec[2]) *
            (PsiX*ch2_vec[0] +PsiY*ch2_vec[1] + PsiZ*ch2_vec[2])) +
            ((ChiX*ch1_vec[0] +ChiY*ch1_vec[1] + ChiZ*ch1_vec[2]) *
            (ChiX*ch2_vec[0] +ChiY*ch2_vec[1] + ChiZ*ch2_vec[2]))) *
            np.exp(2*np.pi*1j*f*(OmgX*x_vec[0] + OmgY*x_vec[1] +
            OmgZ*x_vec[0])/vs) * (dtheta*dphi*3/(8*np.pi))))
    return gammas, ff

def orf_p_directional(ch1_vec, ch2_vec, det1_loc, det2_loc, vp, f, thetamesh=1,
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
    gammas = (np.sin(THETAS) * (OmgX*ch1_vec[0] +
            OmgY*ch1_vec[1] + OmgZ*ch1_vec[2]) * (OmgX*ch2_vec[0] +
            OmgY*ch2_vec[1] + OmgZ*ch2_vec[2]) *
            np.exp(2*np.pi*1j*f*(OmgX*x_vec[0] + OmgY*x_vec[1] +
                OmgZ*x_vec[0])/vp) * (3/(4*np.pi)))
    return gammas, phis, thetas

def orf_s_directional(ch1_vec, ch2_vec, det1_loc, det2_loc, vs, f, thetamesh=1,phimesh=1):
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
    gammas = (np.sin(THETAS) *
        (((PsiX*ch1_vec[0] +PsiY*ch1_vec[1] + PsiZ*ch1_vec[2]) *
        (PsiX*ch2_vec[0] +PsiY*ch2_vec[1] + PsiZ*ch2_vec[2])) +
        ((ChiX*ch1_vec[0] +ChiY*ch1_vec[1] + ChiZ*ch1_vec[2]) *
        (ChiX*ch2_vec[0] +ChiY*ch2_vec[1] + ChiZ*ch2_vec[2]))) *
        np.exp(2*np.pi*1j*f*(OmgX*x_vec[0] + OmgY*x_vec[1] +
        OmgZ*x_vec[0])/vs) * (3/(8*np.pi)))
    return gammas,phis,thetas

def ccStatReadout_s_wave(Y, sigma, ch1_vec, ch2_vec, det1_loc, det2_loc, vs, thetamesh=1,phimesh=1):
    """
    Read out sky map for S-wave based on cross-correlation spectrum.
    Takes an inverse FFT of the frequency series to get something
    in terms of time-delays and then reads off
    the directions based on the time delays.
    Parameters
    ----------
    Y : `array-like`
        cross-correlation spectrum
    sigma : `array-like`
        sigma on Y
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
    print THETAS[:2,:2], THETAS.shape
    print PHIS[:2,:2], PHIS.shape
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
    gamma0 = (np.sin(THETAS) *
        (((PsiX*ch1_vec[0] +PsiY*ch1_vec[1] + PsiZ*ch1_vec[2]) *
        (PsiX*ch2_vec[0] +PsiY*ch2_vec[1] + PsiZ*ch2_vec[2])) +
        ((ChiX*ch1_vec[0] +ChiY*ch1_vec[1] + ChiZ*ch1_vec[2]) *
        (ChiX*ch2_vec[0] +ChiY*ch2_vec[1] + ChiZ*ch2_vec[2]))))
    tau = np.abs(((OmgX*x_vec[0] + OmgY*x_vec[1] +
        OmgZ*x_vec[0])/vs))
    # inverse fft
    Y_tau = Y.ifft()
    readout = interp1d(Y_tau.times, Y_tau.value)
    vals = readout(tau)
    ccStat = (vals*gamma0)
    return ccStat,sigma**2,phis,thetas
def ccStatReadout_p_wave(Y, sigma, ch1_vec, ch2_vec, det1_loc, det2_loc, vs, thetamesh=1,phimesh=1):
    """
    Read out sky map for S-wave based on cross-correlation spectrum.
    Takes an inverse FFT of the frequency series to get something
    in terms of time-delays and then reads off
    the directions based on the time delays.
    Parameters
    ----------
    Y : `array-like`
        cross-correlation spectrum
    sigma : `array-like`
        sigma on Y
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
    print THETAS[:2,:2], THETAS.shape
    print PHIS[:2,:2], PHIS.shape
    dtheta = thetas[1] - thetas[0]
    dphi = phis[1] - phis[0]
    # much faster to vectorize things
    OmgX = np.sin(THETAS)*np.cos(PHIS)
    OmgY = np.sin(THETAS)*np.sin(PHIS)
    OmgZ = (np.cos(THETAS))
    # get vector perp to Omega
    # take cross product of omega and psi to
    # get third vector
    gamma0 = (np.sin(THETAS) * (OmgX*ch1_vec[0] +
            OmgY*ch1_vec[1] + OmgZ*ch1_vec[2]) * (OmgX*ch2_vec[0] +
            OmgY*ch2_vec[1] + OmgZ*ch2_vec[2]))
    tau = np.abs(((OmgX*x_vec[0] + OmgY*x_vec[1] +
        OmgZ*x_vec[0])/vs))
    # inverse fft
    Y_tau = Y.ifft()
    readout = interp1d(Y_tau.times, Y_tau.value)
    vals = readout(tau)
    ccStat = (vals*gamma0)
    return ccStat,sigma**2,phis,thetas

def orf_p_sph(l,m,ch1_vec, ch2_vec, det1_loc, det2_loc, vp, ff=None, thetamesh=1,
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
    Y_LM = sph_harm(m, l, PHIS, THETAS)
    if ff is None:
        ff = np.logspace(-2,1,num=100)
    gammas = np.zeros(ff.size)
    for ii,f in enumerate(ff):
        gammas[ii] = np.sum(np.sum(Y_LM * np.sin(THETAS) * (OmgX*ch1_vec[0] +
            OmgY*ch1_vec[1] + OmgZ*ch1_vec[2]) * (OmgX*ch2_vec[0] +
            OmgY*ch2_vec[1] + OmgZ*ch2_vec[2]) *
            np.exp(2*np.pi*1j*f*(OmgX*x_vec[0] + OmgY*x_vec[1] +
                OmgZ*x_vec[0])/vp) * (dtheta*dphi*3/(4*np.pi))))
    return gammas, ff

