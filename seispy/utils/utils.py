import numpy as np

def set_channel_vector(letter):
    if letter.upper() == 'Z':
        return [0, 0, 1]
    if letter.upper() == 'N':
        return [0, 1, 0]
    if letter.upper() == 'E':
        return [1, 0 ,0]
def _calc_travel_time(delta_vec, Omega, v):
    """
    calculate travel time between two stations assuming
    the ray takes a straight line.
    """
    delta_vec = np.array(delta_vec)
    Omega = np.array(Omega)
    dt = np.dot(delta_vec, Omega)/v
    return dt

def calc_travel_time(delta_vec, OMEGA, v):
    """
    Returns light travel time for many directions on the sky.

    Parameters
    ----------
    delta_vec : `numpy.ndarray`, list
        vector for detector separation. n x 1 matrix.
    OMEGA : `numpy.ndarray`
        matrix of sky direction unit vectors. m x n matrix.
        m is number of sky directions, n is dimension of space.
    v : `float`
        velocity of seismic wave

    Returns
    -------
    dt : `numpy.ndarray`
        list of travel times for each sky direction
    """
    delta_vec = np.array(delta_vec)
    OMEGA = np.array(OMEGA)
    omg_shape = OMEGA.shape
    if not(delta_vec.size == omg_shape[1]):
        raise ValueError('Delta vec should have same dimension as columns of\
        Omega vec')
    dt = np.array([_calc_travel_time(delta_vec, OMEGA[ii],v) for ii in
        range(omg_shape[0])])
    return dt
