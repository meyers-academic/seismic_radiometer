from __future__ import division
from matplotlib import use
use("agg")
from gwpy.plotter.core import Plot
import numpy as np
from matplotlib.projections import register_projection
from matplotlib.projections import AitoffAxes
from plot import FixedMollweideAxes

try:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
    from mpl_toolkits.axes_grid import make_axes_locatable

class RecoveryMapAxes(AitoffAxes):
    """
    axes for recovery matrix"""
    def __init__(self, *args, **kwargs):
        super(RecoveryMapAxes, self).__init__(*args, **kwargs)

    def add_map(self, rec_map, **kwargs):
        dtheta = rec_map.thetas[1] - rec_map.thetas[0]
        dphi = rec_map.phis[1] - rec_map.phis[0]
        # shift from 0,2pi to -pi,pi
        phis = rec_map.phis.copy()
        newphis = np.zeros(phis.size)
        newphis[phis>np.pi] = phis[phis > np.pi] - 2 * np.pi
        newphis[phis<=np.pi] = phis[phis<=np.pi]
        args = np.argsort(newphis)
        print newphis[args]
        self.pcolormesh(newphis[args] - dphi/2, np.pi / 2 - rec_map.thetas +
                        dtheta/2, rec_map.data[args,:].T,
                        **kwargs)

register_projection(RecoveryMapAxes)

class RecoveryMapPlot(Plot):
    """
    recovery map plot
    """
    _DefaultAxesClass = RecoveryMapAxes
    def __init__(self, recovery_map, contour=True, cmap='viridis',vmin=None,
            vmax=None, conf=None, **kwargs):
        super(RecoveryMapPlot, self).__init__(**kwargs)
        dtheta = recovery_map.thetas[1] - recovery_map.thetas[0]
        dphi = recovery_map.phis[1] - recovery_map.phis[0]
        if conf is None:
            conf = 0.5
        #conf_map, phil, thetal = recovery_map.get_contour(conf)
        self.recovery_map = recovery_map
        ax = self.gca()
        ax.add_map(recovery_map, **{'cmap':cmap,'vmin':vmin,'vmax':vmax})
        ax.grid()
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        ax.tick_params(axis='x',colors='white')
        cbar = self.add_colorbar(label=r'amplitude [$\textrm{m}^2$]')
        cbar.ax.tick_params(labelsize=10)
        self.cbar=cbar
        S_sorted = np.sort(recovery_map.data.flatten())[::-1]
        S_cdf = S_sorted.cumsum() / S_sorted.cumsum()[-1]
        #print S_sorted[S_cdf < conf]
        contour_val = np.min(np.real(S_sorted[S_cdf < conf]))
        total_power = np.sum(np.real(S_sorted[S_cdf < conf]))
        ax.set_title('Total power in spot is: %4.2e m$^2$ / Hz' % total_power,
                y=1.08)
        if contour:
            phis = recovery_map.phis.copy()
            newphis = np.zeros(phis.size)
            newphis[phis>np.pi] = phis[phis > np.pi] - 2 * np.pi
            newphis[phis<=np.pi] = phis[phis<=np.pi]
            ax.contour(newphis - dphi/2, np.pi / 2 - recovery_map.thetas + dtheta / 2,
                   recovery_map.data.T, colors='k', linewidth=4,
                   levels=[contour_val])
