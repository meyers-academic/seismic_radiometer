from __future__ import divide
import numpy as np
from ..trace import Trace

def simulate_pwave(A, phi, theta, frequency, time, Fs=100):
    """
    simulate s-wave in a certain direction

    Parameters
    ----------
    direction : TODO
    frequency : TODO
    time : TODO
    Fs : TODO

    Returns
    -------
    E : TODO
    N : TODO
    Z : TODO
    """
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    src_dir = np.array([cphi*ctheta, sphi*ctheta, stheta])
    times = np.arange(0, time, 1/Fs)
    amp = A * np.cos(2*np.pi*frequency*times)
    return src_dir[0]*amp, src_dir[1]*amp, src_dir[2]*amp
