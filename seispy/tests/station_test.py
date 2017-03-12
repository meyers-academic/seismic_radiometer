from __future__ import division
import matplotlib

matplotlib.use('agg')
import unittest
from ..station import *
import numpy.testing as npt


class TestStation(unittest.TestCase):
    def test_spiral(self):
        array = spiral(5, radius=1, height=0, n_per_turn=4)
        npt.assert_almost_equal(array[0], [1, 0, 0])
        npt.assert_almost_equal(array[4], [1, 0, 0])

    def test_homestake(self):
        # check a few cases
        arr = homestake('underground')
        # noinspection PyTypeChecker,PyTypeChecker
        npt.assert_almost_equal(len(arr.keys()), 15)
        arr = homestake('surface')
        # noinspection PyTypeChecker,PyTypeChecker
        npt.assert_almost_equal(len(arr.keys()), 9)
        arr = homestake('all')
        # noinspection PyTypeChecker,PyTypeChecker
        npt.assert_almost_equal(len(arr.keys()), 24)
        arr = homestake()
        # noinspection PyTypeChecker,PyTypeChecker
        npt.assert_almost_equal(len(arr.keys()), 24)
        arr = homestake(['DEAD', 'YATES'])
        # noinspection PyTypeChecker,PyTypeChecker
        npt.assert_almost_equal(len(arr.keys()), 2)


if __name__ == "__main__":
    unittest.main()
