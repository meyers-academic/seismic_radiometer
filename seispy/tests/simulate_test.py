from __future__ import division
import matplotlib
matplotlib.use('agg')
import unittest
import matplotlib.pyplot as plt
from ..simulate import *
import numpy.testing as npt
from gwpy.frequencyseries import FrequencySeries
import numpy as np
STATIONS = {0:[0,0,0]}
A = 10
PHI = 0
THETA = np.pi / 2
PSI = 0
FF = 1
DURATION = 20
SAMPLE_FREQ = 100
VEL = 3000
NOISE = 0
PHASE = 0

class TestSimulations(unittest.TestCase):
    def test_p_wave(self):
        # should be only data in E direction when we prop. in E direction.
        data = simulate_pwave(STATIONS, A, PHI, THETA, FF, DURATION,
                Fs=SAMPLE_FREQ, c=VEL, noise_amp=NOISE, phase=PHASE)
        npt.assert_array_almost_equal(data[0]['E'].value,
                A*np.sin(2*np.pi*FF*np.arange(0,DURATION*SAMPLE_FREQ)/SAMPLE_FREQ))
        npt.assert_array_almost_equal(data[0]['N'].value,
                np.zeros(SAMPLE_FREQ*DURATION))
        npt.assert_array_almost_equal(data[0]['Z'].value,
                np.zeros(SAMPLE_FREQ*DURATION))

    def test_s_wave(self):
        # should only be in N direction if we prop. in E direction.
        data = simulate_swave(STATIONS, A, PHI, THETA, PSI, FF, DURATION,
                Fs=SAMPLE_FREQ, c=VEL, noise_amp=NOISE, phase=PHASE)
        npt.assert_array_almost_equal(data[0]['N'].value,
                A*np.sin(2*np.pi*FF*np.arange(0,DURATION*SAMPLE_FREQ)/SAMPLE_FREQ))
        npt.assert_array_almost_equal(data[0]['E'].value,
                np.zeros(SAMPLE_FREQ*DURATION))
        npt.assert_array_almost_equal(data[0]['Z'].value,
                np.zeros(SAMPLE_FREQ*DURATION))

    def test_get_pol_directions(self):
        # wave travels in E direction
        # psi=0 => polarization should be in N.
        dx, dy, dz = get_polarization_coeffs(0, np.pi / 2, 0)
        npt.assert_almost_equal(dx, 0)
        npt.assert_almost_equal(dy, 1)
        npt.assert_almost_equal(dz, 0)
        # wave travels in N direction
        # just shift your right hand around. it'll make sense.
        # psi=0 => polarization should be in -E.
        dx, dy, dz = get_polarization_coeffs(np.pi/2, np.pi / 2, 0)
        npt.assert_almost_equal(dx, -1)
        npt.assert_almost_equal(dy, 0)
        npt.assert_almost_equal(dz, 0)
        # wave travels in Z direction
        # psi=0 => polarization should be in E.
        dx, dy, dz = get_polarization_coeffs(0, 0, 0)
        npt.assert_almost_equal(dx, 0)
        npt.assert_almost_equal(dy, 1)
        npt.assert_almost_equal(dz, 0)
if __name__=="__main__":
    unittest.main()

