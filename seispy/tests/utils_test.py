import matplotlib
matplotlib.use('agg')
import unittest
from ..utils import *
import matplotlib.pyplot as plt
import numpy.testing as npt
from gwpy.frequencyseries import FrequencySeries

class TestOrfs(unittest.TestCase):
    """
    test overlap reduction functions for s and p-waves
    """
    def test_orfp(self):
        # get them for colocated detectors
        # should be all ones
        gamma, freqs =\
            orf_p(np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),1000)
        test = np.ones(gamma.size)
        npt.assert_array_almost_equal(gamma, test)

    def test_orfs(self):
        # get them for colocated detectors
        # should be all ones
        gamma, freqs =\
            orf_s(np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),1000)
        test = np.ones(gamma.size)
        npt.assert_array_almost_equal(gamma, test, decimal=4)

    def test_orfp_dir(self):
        gamma, phis, thetas=\
            orf_p_directional(np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),1000,
                    1)

    def test_orfs_dir(self):
        gamma, phis, thetas=\
            orf_s_directional(np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),1000,
                    1)

    def test_orf_p_sph(self):
        # get them for colocated detectors
        # should be all ones
        gamma, freqs =\
            orf_p_sph(0,0,np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1000,0,0]),1000,
                    ff=np.logspace(-2,1,1000))

    def test_ccStatReadout_s_wave(self):
        Y = FrequencySeries(np.random.randn(1000)+1j*np.random.randn(1000),
                frequencies=np.arange(1000)*0.01)
        ccStat,sigma, phis, thetas=\
            ccStatReadout_s_wave(Y, Y,
                    np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1,0,0]),np.asarray([1000,200,312]),1000)

if __name__=="__main__":
    unittest.main()

