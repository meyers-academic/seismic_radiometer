from __future__ import division
from matplotlib import use
use("agg")
from gwpy.plotter.core import Plot
import numpy as np
from matplotlib.projections import register_projection
from matplotlib.projections import AitoffAxes

try:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
    from mpl_toolkits.axes_grid import make_axes_locatable


class RecoveryMapAxes(AitoffAxes):
    """
    axes for recovery matrix"""
    def __init__(self, *args, **kwargs):
        super(RecoveryMapAxes, self).__init__(*args, **kwargs)

    def add_map(self, rec_map):
        dtheta = rec_map.thetas[1] - rec_map.thetas[0]
        dphi = rec_map.phis[1] - rec_map.phis[0]
        self.pcolormesh(rec_map.phis - np.pi - dphi/2, np.pi / 2 - rec_map.thetas + dtheta/2, rec_map.data.T,
        cmap='viridis')

register_projection(RecoveryMapAxes)

class RecoveryMapPlot(Plot):
    """
    recovery map plot
    """
    _DefaultAxesClass = RecoveryMapAxes
    def __init__(self, recovery_map, contour=True, **kwargs):
        super(RecoveryMapPlot, self).__init__(**kwargs)
        dtheta = recovery_map.thetas[1] - recovery_map.thetas[0]
        dphi = recovery_map.phis[1] - recovery_map.phis[0]
        try:
            conf = kwargs.pop('conf')
        except KeyError:
            conf = 0.5
        conf_map, phil, thetal = recovery_map.get_contour(conf)
        self.recovery_map = recovery_map
        ax = self.gca()
        ax.add_map(recovery_map)
        ax.grid()
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        ax.tick_params(axis='x',colors='white')
        cbar = self.add_colorbar(label=r'amplitude [$\textrm{m}^2$]')
        cbar.ax.tick_params(labelsize=10)
        if contour:
            ax.contour(recovery_map.phis - np.pi - dphi/2, np.pi / 2 - recovery_map.thetas + dtheta / 2,
                   conf_map.T, colors='k', linewidth=4, levels=[0])

