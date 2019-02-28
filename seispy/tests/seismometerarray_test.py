from __future__ import division
import matplotlib

matplotlib.use('agg')
import unittest
from ..station import SeismometerArray
import numpy.testing as npt
import numpy as np

STATIONS = {0: [0, 0, 0]}
A = 10
PHI = 0
THETA = np.pi / 2
PSI = 0
FF = 1
ALPHA = 1000
EPSILON = 0.1
DURATION = 20
SAMPLE_FREQ = 100
VEL = 3000
NOISE = 0
PHASE = 0


# noinspection PyTypeChecker
class TestSeismometerArray(unittest.TestCase):
    def test_p_wave(self):
        """
        test p-wave for simple case of parameters defined at top of this file
        """
        data = SeismometerArray.initialize_all_good(STATIONS, DURATION, chans_type='fast_chans', start_time=0)
        data.add_p_wave(A, PHI, THETA, FF, DURATION)
        print data[0]['HHE'].psd().sum()
        # should be only data in E direction when we prop. in E direction.
        npt.assert_array_almost_equal(data[0]['HHE'].value,
                                      A * np.sin(2 * np.pi * FF * np.arange(0, DURATION * SAMPLE_FREQ) / SAMPLE_FREQ))
        npt.assert_array_almost_equal(data[0]['HHN'].value,
                                      np.zeros(SAMPLE_FREQ * DURATION))
        npt.assert_array_almost_equal(data[0]['HHZ'].value,
                                      np.zeros(SAMPLE_FREQ * DURATION))

    def test_s_wave(self):
        """
        test s-wave for simple case of parameters defined at top of this file.
        """
        # should only be in N direction if we prop. in E direction.
        data = SeismometerArray.initialize_all_good(STATIONS, DURATION, chans_type='fast_chans', start_time=0)
        data.add_s_wave(A, PHI, THETA, PSI, FF, DURATION)
        # noinspection PyTypeChecker
        npt.assert_array_almost_equal(data[0]['HHN'].value,
                                      A * np.sin(2 * np.pi * FF * np.arange(0, DURATION * SAMPLE_FREQ) / SAMPLE_FREQ))
        npt.assert_array_almost_equal(data[0]['HHE'].value,
                                      np.zeros(SAMPLE_FREQ * DURATION))
        npt.assert_array_almost_equal(data[0]['HHZ'].value,
                                      np.zeros(SAMPLE_FREQ * DURATION))

#     def test_r_wave(self):
#         """
#         test r-wave for simple case of parameters
#         defined at top of this file
#         """
#         # should only be in N direction if we prop. in E direction.
#         data = SeismometerArray.initialize_all_good(STATIONS, DURATION, chans_type='fast_chans', start_time=0)
#         data.add_r_wave(A, PHI, THETA, EPSILON, ALPHA, FF, DURATION)
#         # noinspection PyTypeChecker
#         npt.assert_array_almost_equal(data[0]['HHE'].value,
#                                       A * np.cos(2 * np.pi * FF * np.arange(0, DURATION * SAMPLE_FREQ) / SAMPLE_FREQ))
#         # don't forget to multiply by epsilon!
#         # noinspection PyTypeChecker
#         npt.assert_array_almost_equal(data[0]['HHZ'].value,
#                                       -EPSILON * A * np.sin(
#                                           2 * np.pi * FF * np.arange(0, DURATION * SAMPLE_FREQ) / SAMPLE_FREQ))
#         npt.assert_array_almost_equal(data[0]['HHN'].value,
#                                       np.zeros(SAMPLE_FREQ * DURATION))
# 

if __name__ == "__main__":
    unittest.main()
