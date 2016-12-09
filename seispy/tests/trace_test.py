from ..trace import Trace, fetch
import unittest

EPOCH_START = 1112720556
EPOCH_END = 1112720756
BASE_DIR = 'data/'
CHANNEL = 'D4850:HHZ'
EXPECTED_SRATE = '100'

class TraceTest(unittest.TestCase):
    def test_read_frame():
        data = fetch(EPOCH_START, EPOCH_END, CHANNEL, framedir=BASE_DIR)
        self.assertTrue(isinstance(data, Trace))
        self.assertTrue(data.times.value[0] == EPOCH_START)
        self.assertTrue(data.times.value[-1] == EPOCH_END)
        self.assertTrue(data.dx.value == 1./EXPECTED_SRATE)

if __name__=="__main__":
    unittest.main()
