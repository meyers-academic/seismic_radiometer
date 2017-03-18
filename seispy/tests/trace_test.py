from ..trace import Trace, fetch
import unittest
import numpy as np
import numpy.testing as npt
import obspy


EPOCH_START = 1125384593
EPOCH_END = 1125384693
BASE_DIR = 'seispy/tests/data/'
CHANNEL = 'D4850:HHZ'
EXPECTED_SRATE = 100



class TraceTest(unittest.TestCase):
    def test_read_frame(self):
        data = fetch(EPOCH_START, EPOCH_END, CHANNEL, framedir=BASE_DIR)
        self.assertTrue(isinstance(data, Trace))
        self.assertTrue(data.times.value[0] == EPOCH_START)
        self.assertTrue(data.times.value[-1] == EPOCH_END - data.dx.value)
        self.assertTrue(data.dx.value == 1. / EXPECTED_SRATE)

    def test_vel2disp(self):
        """
        Check velocity to displacement method

        Returns
        -------

        """
        # add a sine-wave of frequency 2
        data = Trace(np.cos(2 * np.pi * 2 * np.arange(0, 100, 0.01)), sample_rate=100)
        # generate an ASD and get its max
        asd = data.asd()
        amax = np.argmax(asd)
        # change from velocity to displacement.
        # means we divide by frequency
        new_trace = data.vel2disp()
        # get a new ASD
        asd2 = new_trace.asd()
        # check that we've gone down by factor of 2 at 2 Hz
        npt.assert_almost_equal(0.5, asd2[amax].value/asd[amax].value)

    def test_gaussian_filter(self):
        """
        test application of a gaussian filter

        Returns
        -------

        """
        # add a sine-wave
        data = Trace(np.cos(2*np.pi * 2 * np.arange(0, 100, 0.01)), sample_rate=100)
        # add in some noise so we can check our plots
        data += np.random.randn(data.size)
        # take asd and its maximum
        asd_original = data.asd()
        amax = np.argmax(asd_original)
        # apply a gaussian filter
        newtrace = data.gaussian_filter(2, 1)
        # take asd after filtering
        asd_filtered = newtrace.asd()
        # make sure peak power is unchanged
        npt.assert_almost_equal(asd_original[amax].value,
                                asd_filtered[amax].value,
                                decimal=5)

    def test_toobspy(self):
        # add a sine-wave
        data = Trace(np.cos(2*np.pi * 2 * np.arange(0, 100, 0.01)),
                     sample_rate=100, channel=CHANNEL)
        data2 = data.to_obspy()
        self.assertTrue(isinstance(data2,
                                   obspy.core.trace.Trace))

if __name__ == "__main__":
    unittest.main()
