from seispy.trace import (Trace, fetch)

EPOCH_START = 1112720556
EPOCH_END = 1112720756
BASE_DIR = 'data/'
CHANNEL = 'D4850:HHZ'
EXPECTED_SRATE = '100'

data = fetch(EPOCH_START, EPOCH_END, CHANNEL, framedir=BASE_DIR)
if not isinstance(data, Trace):
	raise TypeError('read_frame should output Trace class')
if not data.times.value[0] == EPOCH_START:
	raise ValueError('start of data read out should be %d but is actually %d' % (EPOCH_START, data.times.value[0]))
if not data.times.value[-1] == EPOCH_END:
	raise ValueError('start of data read out should be %d but is actually %d' % (EPOCH_END, data.times.value[-1]))
if not data.dx.value == 1./EXPECTED_SRATE:
	raise ValueError('data sample rate should be %4.2f but is actually %4.2f' % (1./EXPECTED_SRATE, data.dx.value))

print 'TEST PASSED!'