from station import Station
from collections import OrderedDict


class StationArray(OrderedDict):
    """
    StationArray class

    >>> frame = 'M-SEISMIC-1125384593-4096.gwf'
    >>> station_names = ['DEAD','D4850']
    >>> arr = StationArray(station_names, st, et)
    >>> print arr
    >>> arr[station_name[0]]['Z']
    """

    def __init__(self, st, et, stations):
        super(StationArray, self).__init__()
        self.st = st
        self.et = et
        self.stations = stations
        for station in stations:
            self[station] = Station(st, et, station)

    def coherence(self, fftlength=None, window='hanning',
                  **kwargs):
        COH = OrderedDict()
        for ii in range(len(self.keys())):
            for jj in range(ii, len(self.keys())):
                if ii == jj:
                    continue
                newkey = self.keys()[ii] + '-' + self.keys()[jj]
                key1 = self.keys()[ii]
                key2 = self.keys()[jj]
                COH[newkey] = self[key1]['Z'].coherence(self[key2]['Z'],
                                                        fftlength=fftlength,
                                                        window=window,
                                                        stacktype='ts',
                                                        **kwargs)
        return COH
