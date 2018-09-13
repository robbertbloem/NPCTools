"""
Test change 2
"""


from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import imp

import scipy
import numpy
import matplotlib 
import matplotlib.pyplot as plt

import struct
import re

import NPCTools.Resources.CommonFunctions as CF
import NPCTools.Resources.Mathematics as MATH
import NPCTools.Resources.Equations as EQ
import NPCTools.Resources.StreakFunctions as SF


imp.reload(CF)
imp.reload(SF)


class StreakTools():
    
    def __init__(self, path, filename, debug = 0):
        self.path = path
        self.filename = filename
        self.paf = path + filename
        self.debug = debug
    
    
    def import_streak(self):
        self.z, self.w_axis, self.t_axis, self.meta = SF.import_streak(self.paf)


    def print_meta(self): 
        for k, v in self.meta.items():
            print("%s" % k.upper())
            for x, y in v.items():
                print("\t%25s = %s" % (x, y))
    
    def get_w_axis(self):
        return self.w_axis
    
    def get_w_axis_label(self):
        return "Wavelength (nm)"

    def get_time_axis(self):
        return self.t_axis
        
    def get_time_axis_label(self):
        s = ""
        return s


    def get_area(self, w_range = [0,-1], t_range = [0,-1], frame = "data", store = False):

        _w_axis, _t_axis, _z = SF.get_area(self.w_axis, self.t_axis, self.z, w_range = w_range, t_range = t_range, frame = frame, debug = self.debug)
        
        if store:
            self.w_axis = _w_axis
            self.t_axis = _t_axis
            self.z = _z
        
        return _w_axis, _t_axis, _z


    def get_average(self, axis = "w", range = [0,-1], frame = "data", return_cumulative = False):
        """
        INPUTS:
        w_axis: wavelength axis
        t_axis: time axis
        z: data
        axis: 'w' or 't', the axis you want to keep. Average w means that the t-axis is averaged. 
        range: start and end of the range to be averaged.
        frame: 'data' or 'index'
        """
        
        _z = SF.get_average(self.w_axis, self.t_axis, self.z, axis = axis, range = range, frame = frame, return_cumulative = return_cumulative, debug = self.debug)
        
        return _z






    def plot_streak(self, ax = False, **kwargs):

        """
        x_range, y_range (list with 2 elements, default [None, None]): 
        ax: axis object, if False, a new figure will be created.
    
        kwargs:
        - z_range: [min, max]
    
        """

        if ax == False:
            fig = plt.figure()
            ax = fig.add_subplot(111)
    
#         self.z = -self.z
    
        if "z_range" in kwargs:
            levels = numpy.linspace(kwargs["z_range"][0], kwargs["z_range"][1])
            kwargs["levels"] = levels
        
        ax.contourf(self.w_axis, self.t_axis, self.z, **kwargs)
    
        ax.set_xlabel("Wavelength (%s)" % self.meta["Scaling"]["ScalingXUnit"])
        ax.set_ylabel("Time (%s)" % self.meta["Scaling"]["ScalingYUnit"])

        if "x_range" in kwargs:
            ax.set_xlim(kwargs["x_range"])

        if "y_range" in kwargs:
            ax.set_ylim(kwargs["y_range"])

        return ax


    def fit_lifetime(self, w_range = [0,-1], t_range = [0,-1], frame = "data", flag_weigh_laser = False, fit_type = "single_exp", exp_hints = -1):

        _t_axis, _z, y_fit, A_out = SF.fit_lifetime(self.w_axis, self.t_axis, self.z, w_range = w_range, t_range = t_range, frame = frame, flag_weigh_laser = flag_weigh_laser, fit_type = fit_type, debug = self.debug, exp_hints = exp_hints)
        
        return _t_axis, _z, y_fit, A_out
        
        
    def make_lifetime_plot(self, _t_axis, _z, y_fit, A_out, fig = False, ax_decay = False, ax_residue = False):
        
        if fig == False:
            figsize = (8.0, 6.0)
            coords = [(0.0875, 0.4167, 0.875, 0.5), (0.0875, 0.08333, 0.875, 0.3333)]
            n_fig = 1
            n_ax = len(coords)
            fig = [0] * n_fig
            ax = [0] * n_fig
            for fig_i in range(n_fig):    
                fig[fig_i] = plt.figure(figsize = figsize)
                ax[fig_i] = [0] * n_ax
                for ax_i in range(n_ax):
                    ax[fig_i][ax_i] = fig[fig_i].add_axes(coords[ax_i])
                ax_decay = ax[0][0]
                ax_residue = ax[0][1]
                
        if ax_decay:
            ax_decay.plot(_t_axis, _z)
            ax_decay.plot(_t_axis, y_fit)
        
        if ax_residue:
            ax_residue.plot(_t_axis, _z - y_fit)
                
            
        



if __name__ == "__main__": 
    pass