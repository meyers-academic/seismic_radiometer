from __future__ import division
from matplotlib import use,rc
use('agg')
rc('text',usetex=True)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from collections import OrderedDict
import numpy as np

class StationArray(OrderedDict):
    """
    Station array class. Inherits from dict"""
    def __init__(self):
        super(StationArray, self).__init__()

    def shift_origin(self,origin=None):
        """
        shift the origin of the StationArray (changes StationArray in place)

        Parameters
        ----------
        origin: name of station to use as origin.
           if no origin is specified, this will use the first
           origin in the list.
        """
        # if there are no elements, just return
        if len(self)==0:
            return
        # if origin is none, get first element of array and call that the origin
        if origin==None:
            origin=self.keys()[0]
            origin_x,origin_y,origin_z=self[origin][0],self[origin][1],self[origin][2]
            origin_coords=[origin_x,origin_y,origin_z]
        # otherwise, get the coordinates of the user specified origin
        else:
            try:
                origin_x,origin_y,origin_z=self[origin][0],self[origin][1],self[origin][2]
                origin_coords=[origin_x,origin_y,origin_z]
            except:
                print 'Unknown station used as origin'
                raise
        # shift coordinates of all array elements to the new origin
        for station in self.keys():
            for ii in range(3):
                self[station][ii]=self[station][ii]-origin_coords[ii]
                

    def plot(self):
        """
        plot station array

        Parameters
        ----------
        *args : TODO
        **kwargs : TODO

        Returns
        -------
        TODO

        """
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        for key in stations.keys():
            ax.scatter(stations[key][0],stations[key][1],stations[key][2])
        ax.set_zlabel('Depth')
        ax.set_ylabel('North')
        ax.set_xlabel('East')
        return plt


def spiral(N, radius=1000, height=100, n_per_turn=10, offset=0, origin=None):
    """
    return a spiral array pattern

    Parameters
    ----------
    N : `int`
        number of stations
    radius : `float`, optional
        radius of spiral circles
    height : `float`, optional
        height step from one station to the next.
    n_per_turn : `float`, optional
        number of stations in a circle (when projected onto z-axis)
    origin: `int` optional
        station to use as an origin. defaults to using first station in list.

    Returns
    -------
    stations : :class:`seispy.station.StationArray`
        an ordered dict of stations
    """
    stations = StationArray()
    for ii in range(N):
        stations[ii] = offset + radius * np.array([np.cos((1/n_per_turn)*2*np.pi*ii),
            np.sin((1/n_per_turn)*2*np.pi*ii), ii*height/radius])

    # shift origin
    stations.shift_origin(origin)

    return stations

def homestake(keylist=None,origin=None):
    """
    Return a dict with homestake array

    Parameters
    ----------
    stationlist : `str`, optional
        list of stations to return in StationArray object. Default
        is to return all stations. Will also return all if 'All' is given.
        Will return surface stations of 'surface' is given.
        Will return underground stations if 'underground' is given.
    origin : `str`, optional
        station which will be used as origin (None defaults to using first element in array)
    """
    stations = StationArray()
    xyz_list = {'DEAD': [599316.496208, 4915135.795515, 1498],
                'LHS': [597654.698101, 4911177.756239, 1684],
                'ORO': [599454.495625, 4910782.727115, 1543],
                'ROSS': [599029.335575, 4910954.029719, 1628],
                'RRDG': [598383.512647, 4912544.118885, 1677],
                'SHL': [602888.695761, 4907880.584008, 1772],
                'TPK': [595816.121345, 4910428.385396, 1740],
                'WTP': [600233.873095, 4911950.092981, 1555],
                'YATES': [599503.542430, 4911750.052702, 1625],
                '300': [599082.941268, 4911099.273511, 1505.1],
                '800': [598921.845912, 4911207.932696, 1350],
                '1700': [598954.462337, 4911686.159936, 1073.8],
                'A2000': [598644.641050, 4911614.812799, 983.3],
                'B2000': [598791.293544, 4911405.938094, 983.3],
                'C2000': [598532.222831, 4911668.666166, 983.1],
                'D2000': [598114.948236, 4911851.254988, 983],
                'E2000': [597894.603780, 4912192.359339, 983.1],
                'A4100': [599356.477888, 4910936.776895, 342.5],
                'C4100': [599568.994798, 4911639.949104, 342.3],
                'D4100': [599549.979057, 4910795.291477, 342.4],
                'A4850': [598646.478984, 4910437.175145, 115.2],
                'B4850': [598987.461619, 4911086.715912, 114.9],
                'C4850': [599409.907522, 4911093.130999, 114.6],
                'D4850': [599581.886292, 4911840.127688, 115.2]}
    surface_list = ['DEAD','LHS','ORO','ROSS','RRDG','SHL','TPK','WTP','YATES']
    ug_list = [] # all others
    for key in xyz_list.keys():
        cont=0
        for t in surface_list:
            if t==key:
                cont=1
        if cont==1:
            continue
        else:
            ug_list.append(key)

    if isinstance(keylist,str):
        if keylist.lower() == 'all':
            keylist=xyz_list.keys()
        elif keylist.lower() == 'surface':
            keylist = surface_list
        elif keylist.lower() == 'underground':
            keylist = ug_list
    elif keylist is None:
        keylist = xyz_list.keys()
    for key in keylist:
        stations[key]=xyz_list[key]

    # shift origin
    stations.shift_origin(origin)

    return stations
