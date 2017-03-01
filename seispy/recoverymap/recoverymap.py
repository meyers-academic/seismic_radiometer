from __future__ import division
from matplotlib import use
from ..plotter.recoverymap import RecoveryMapPlot
use('agg')
import numpy as np
import matplotlib.pyplot as plt

class RecoveryMap(object):
    """
    recovery map"""
    def __init__(self, data, thetas, phis, maptype):
        super(RecoveryMap, self).__init__()
        self.data = data
        self.thetas = thetas
        self.phis = phis
        self.maptype = maptype

    def get_contour(self, conf):
        """
        get confidence contour
        """
        map_shape=self.data.shape
        flat_map = self.data.flatten().squeeze()
        flat_map[flat_map < 0] = 0
        args = np.argsort(flat_map)[::-1]
        cdf = flat_map[args].cumsum()
        cdf = cdf/cdf[-1]
        cdf_map = np.zeros(cdf.size)
        cdf_map[args]=cdf
        conf_args = args[cdf<conf]
        Map_conf_percent = np.zeros(cdf.size)
        Map_conf_percent[conf_args]=1
        map_conf = Map_conf_percent.reshape(map_shape)
        return map_conf

    def power_in_conf(self, conf):
        """
        get power in confidence region
        """
        map_shape=self.data.shape
        flat_map = self.data.flatten().squeeze()
        flat_map[flat_map < 0] = 0
        args = np.argsort(flat_map)[::-1]
        cdf = flat_map[args].cumsum()
        cdf = cdf/cdf[-1]
        cdf_map = np.zeros(cdf.size)
        cdf_map[args]=cdf
        conf_args = args[cdf<conf]
        print 'CONF IS:'  + str(conf)
        Map_conf_percent = np.zeros(cdf.size)
        Map_conf_percent[conf_args]=1
        map_conf = Map_conf_percent.reshape(map_shape)
        return np.sqrt(flat_map[conf_args].sum())

    def plot(self, *args, **kwargs):
        return RecoveryMapPlot(self, **kwargs)
