import matplotlib

matplotlib.use('agg')
import unittest
from ..utils.orfs import *
from ..utils.utils import *
import numpy.testing as npt
import numpy as np


class TestOrfs(unittest.TestCase):
    """
    test overlap reduction functions for s and p-waves
    """

    def test_orfp(self):
        # get them for colocated detectors
        # should be all ones
        gamma, freqs = \
            orf_p(np.asarray([1, 0, 0]), np.asarray([1, 0, 0]), np.asarray([1, 0, 0]), np.asarray([1, 0, 0]), 1000)
        test = np.ones(gamma.size)
        npt.assert_array_almost_equal(gamma, test)

    def test_orfs(self):
        # get them for colocated detectors
        # should be all ones
        gamma, freqs = \
            orf_s(np.asarray([1, 0, 0]), np.asarray([1, 0, 0]), np.asarray([1, 0, 0]), np.asarray([1, 0, 0]), 1000)
        test = np.ones(gamma.size)
        npt.assert_array_almost_equal(gamma, test, decimal=4)

    def test_orfp_dir(self):
        # for same detector, same channel, source along channel direction, we
        # expect orf to be 1.
        gamma, phis, thetas = \
            orf_p_directional(np.asarray([1, 0, 0]), np.asarray([1, 0, 0]), np.asarray([1, 0, 0]),
                              np.asarray([1, 0, 0]), 1000,
                              1, thetas=np.array([np.pi / 2, np.pi / 2]),
                              phis=np.array([0, 0]))
        npt.assert_array_almost_equal(gamma.squeeze(),
                                      np.ones(gamma.size).reshape(gamma.shape))

    def test_orfs_dir(self):
        # check with channel aligned with polarization 1
        gamma1, gamma2, phis, thetas = \
            orf_s_directional(np.asarray([0, 1, 0]), np.asarray([0, 1, 0]), np.asarray([1, 0, 0]),
                              np.asarray([1, 0, 0]), 1000,
                              1, thetas=np.array([np.pi / 2, np.pi / 2]),
                              phis=np.array([0, 0]))
        npt.assert_array_almost_equal(gamma1.squeeze(),
                                      np.ones(gamma1.size).reshape(gamma1.shape))
        npt.assert_array_almost_equal(gamma2.squeeze(),
                                      np.zeros(gamma2.size).reshape(gamma2.shape))
        # check with channel aligned with polarization 2
        gamma1, gamma2, phis, thetas = \
            orf_s_directional(np.asarray([0, 0, 1]), np.asarray([0, 0, 1]), np.asarray([1, 0, 0]),
                              np.asarray([1, 0, 0]), 1000,
                              1, thetas=np.array([np.pi / 2, np.pi / 2]),
                              phis=np.array([0, 0]))
        npt.assert_array_almost_equal(gamma2.squeeze(),
                                      np.ones(gamma2.size).reshape(gamma2.shape))
        npt.assert_array_almost_equal(gamma1.squeeze(),
                                      np.zeros(gamma1.size).reshape(gamma1.shape))

        # def test_orf_p_sph(self):
        #     # get them for colocated detectors
        #     # should be all ones
        #     gamma, freqs =\
        #         orf_p_sph(0,0,np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1000,0,0]),1000,
        #                 ff=np.logspace(-2,1,1000))
        #
        # def test_ccStatReadout_s_wave(self):
        #     Y = FrequencySeries(np.random.randn(1000)+1j*np.random.randn(1000),
        #             frequencies=np.arange(1000)*0.01)
        #     ccStat,sigma, phis, thetas=\
        #         ccStatReadout_s_wave(Y, Y,
        #                 np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1000,200,312]),1000)


class TestUtils(unittest.TestCase):
    def test_light_travel_time(self):
        v = 1000
        delta_vec = [0, 0, 0]
        OMEGA = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        dt = calc_travel_time(delta_vec, OMEGA, v)
        npt.assert_array_almost_equal(dt, np.zeros(dt.size))
        delta_vec = [1000, 0, 0]
        dt = calc_travel_time(delta_vec, OMEGA, v)
        npt.assert_array_almost_equal(dt, np.array([1, 0, 0]))


if __name__ == "__main__":
    unittest.main()
