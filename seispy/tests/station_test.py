from __future__ import division
import matplotlib
matplotlib.use('agg')
import unittest
import matplotlib.pyplot as plt
from ..station import *
import numpy.testing as npt
from gwpy.frequencyseries import FrequencySeries
import numpy as np

class TestStation(unittest.TestCase):
    def test_spiral(self):
        array = spiral(5, radius=1, height=0, n_per_turn=4)
        npt.assert_almost_equal(array[0], [1,0,0])
        npt.assert_almost_equal(array[4], [1,0,0])
        print array
    def test_homestake(self):
        # check a few cases
        arr = homestake('underground')
        npt.assert_almost_equal(len(arr.keys()), 15)
        arr = homestake('surface')
        npt.assert_almost_equal(len(arr.keys()), 9)
        arr = homestake('all')
        npt.assert_almost_equal(len(arr.keys()), 24)
        arr = homestake()
        npt.assert_almost_equal(len(arr.keys()), 24)
        arr = homestake(['DEAD','YATES'])
        npt.assert_almost_equal(len(arr.keys()), 2)

if __name__=="__main__":
    unittest.main()

