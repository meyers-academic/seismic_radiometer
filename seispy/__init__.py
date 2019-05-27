constants = dict(km_SI=1000,
                 G_SI=6.673e-11)
default_rwave_parameters = dict(c2=-0.76, c4=-0.69,
                                Nv=-0.68, a1=0.86,
                                a2=0.63, a3=0.63,
                                a4=0.49,
                                alpha=-0.25)


def default_eigenfunction(drp, freq, vr):
    import numpy as np
    V1 = 1j * (1 - drp['c2'])
    V2 = 1j * drp['c2']
    H1 = drp['Nv'] - drp['c4']
    H2 = drp['c4']
    krhomag = 2 * np.pi * freq / vr
    v1 = krhomag * drp['a1']
    v2 = krhomag * drp['a2']
    h1 = krhomag * drp['a3']
    h2 = krhomag * drp['a4']
    return [H1, H2, V1, V2, h1, h2, v1, v2]
