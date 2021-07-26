import numpy as np
import scipy

def pwave_reflection_coefficients(polar_angle, cp, cs):
    """
    Parameters:
    -----------
    incident_angle : float
        polar incident angle in radians
    cp : float
        p-wave velocity
    cs : float
        s-wave velocity
    Returns:
    --------
    """
    incident_angle = np.pi / 2. - polar_angle
    cia = np.cos(incident_angle)
    sia = np.sin(incident_angle)
    reflected_p_angle = -incident_angle
    array_parameter = sia / cp  # often called "p"
    # p = sin(i) / cp = sin(j) / cs
    # => j = arcsin(p * cs)
    reflected_s_angle = np.arcsin(array_parameter * cs)
    srs = np.sin(reflected_s_angle)
    crs = np.cos(reflected_s_angle)
    PP_num = -(cs**-2 - 2 * array_parameter**2) ** 2 + \
        4 * array_parameter**2 * (cia / cp) * (crs / cs)
    PP_denom = (cs**-2 - 2 * array_parameter**2) ** 2 + \
        4 * array_parameter**2 * (cia / cp) * (crs / cs)
    PP = PP_num / PP_denom
    PS_num = 4 * (cp / cs) * array_parameter * (cia / cp) * \
        (cs**-2 - 2 * array_parameter**2)
    PS_denom = (cs**-2 - 2 * array_parameter**2)**2 + \
        4 * array_parameter**2 * (cia / cp) * (crs / cs)
    PS = PS_num / PS_denom
    return PP, PS, reflected_s_angle
