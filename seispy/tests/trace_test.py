from ..trace import Trace, fetch
import unittest

EPOCH_START = 1125384593
EPOCH_END = 1125384693
BASE_DIR = './seispy/tests/data/'
CHANNEL = 'D4850:HHZ'
EXPECTED_SRATE = 100

class TraceTest(unittest.TestCase):
    def test_read_frame(self):
        data = fetch(EPOCH_START, EPOCH_END, CHANNEL, framedir=BASE_DIR)
        self.assertTrue(isinstance(data, Trace))
        self.assertTrue(data.times.value[0] == EPOCH_START)
        self.assertTrue(data.times.value[-1] == EPOCH_END - data.dx.value)
        self.assertTrue(data.dx.value == 1./EXPECTED_SRATE)

if __name__=="__main__":
    unittest.main()
