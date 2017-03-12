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
        if self.maptype=='r':
            idx = np.where(thetas==np.pi / 2.)
            self.data = self.data[:,idx[0]]
            self.thetas=np.asarray([np.pi / 2.])

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
        PHIS,THETAS = np.meshgrid(self.phis, self.thetas)
        phi_rec = PHIS[flat_map.reshape(map_shape).T == np.max(flat_map)]
        theta_rec = THETAS[flat_map.reshape(map_shape).T == np.max(flat_map)]
        min_phi = np.min(PHIS[map_conf.T>0])
        max_phi = np.max(PHIS[map_conf.T>0])
        min_theta = np.min(THETAS[map_conf.T>0])
        max_theta = np.max(THETAS[map_conf.T>0])
        return map_conf, [min_phi, phi_rec[0], max_phi],[min_theta, theta_rec, max_theta]

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
        Map_conf_percent = np.zeros(cdf.size)
        Map_conf_percent[conf_args]=1
        map_conf = Map_conf_percent.reshape(map_shape)
        return np.sqrt(flat_map[conf_args].sum())

    def plot(self, *args, **kwargs):
        return RecoveryMapPlot(self, **kwargs)
