from collections import OrderedDict
from trace import read_frame


class Station(OrderedDict):
    """
    Station class

    initialize station class by:
    >>> frame = 'M-SEISMIC-1125384593-4096.gwf'
    >>> station_name = 'DEAD'
    >>> sta = Station(station_name, st, et)
    >>> print sta
    """

    def __init__(self, st, et, station):
        super(Station, self).__init__()
        self.station = station
        self.st = st
        self.et = et
        self['Z'] = read_frame('M-SEISMIC-1125384593-4096.gwf',
                               station + ':HHZ', st=st, et=et)
        self['N'] = read_frame('M-SEISMIC-1125384593-4096.gwf',
                               station + ':HHN', st=st, et=et)
        self['E'] = read_frame('M-SEISMIC-1125384593-4096.gwf',
                               station + ':HHE', st=st, et=et)
        self.location = self['Z'].location
