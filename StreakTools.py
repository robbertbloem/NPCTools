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


    def get_area(self, w_range = [0,-1], t_range = [0,-1], frame = "data", store = False):

        _w_axis, _t_axis, _z = SF.get_area(self.w_axis, self.t_axis, self.z, w_range = w_range, t_range = t_range, frame = frame, debug = self.debug)
        
        if store:
            self.w_axis = _w_axis
            self.t_axis = _t_axis
            self.z = _z
        
        return _w_axis, _t_axis, _z


    def get_average(self, axis = "w", range = [0,-1], frame = "data"):
        """
        INPUTS:
        w_axis: wavelength axis
        t_axis: time axis
        z: data
        axis: 'w' or 't', the axis you want to keep. Average w means that the t-axis is averaged. 
        range: start and end of the range to be averaged.
        frame: 'data' or 'index'
        """
        
        _z = SF.get_average(self.w_axis, self.t_axis, self.z, axis = axis, range = range, frame = frame, debug = self.debug)
        
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
        print(fit_type)
        _t_axis, _z, y_fit = SF.fit_lifetime(self.w_axis, self.t_axis, self.z, w_range = w_range, t_range = t_range, frame = frame, flag_weigh_laser = flag_weigh_laser, fit_type = fit_type, debug = self.debug, exp_hints = exp_hints)
        
        return _t_axis, _z, y_fit


# def fit_laser_pulse(x, y, z, xrange = [0,-1], yrange = [0,-1], frame = "data", title = "", ax_fit = False, ax_res = False, flag_weigh_laser = False):
# 
#     _z = get_average(x, y, z, axis = "y", range = xrange, frame = frame)
#     
#     if yrange[0] != 0 and yrange[1] != -1:
#         a, b = CF.find_indices_for_range(y, yrange, frame = frame, round = "maximize")
#         
#         _z = _z[a:b]
#         _y = y[a:b]
#     else:
#         _y = numpy.copy(y)
# 
#     # initial try for laser pulse
#     # SD, mean, offset, scale
#     i = numpy.argmax(_z)
#     A_start = [0.1, _y[i], 0, _z[i]]
# 
#     print("FIT RESULTS LASER %s" % (title))
# 
#     if flag_weigh_laser:
# #         _z += 0.00001
#         
#         weight = numpy.copy(_z)
#         temp = numpy.where(weight == 0)[0]
#         weight[temp] += 1e20
#         weight = numpy.sqrt(weight)
#         A_laser, matcov = scipy.optimize.curve_fit(f = gaussian, xdata = _y, ydata = _z, p0 = A_start, sigma = weight, absolute_sigma = True)
#         print("  Weighed")
#     else:
#         A_laser, matcov = scipy.optimize.curve_fit(f = gaussian, xdata = _y, ydata = _z, p0 = A_start)
#         weight = numpy.ones(len(_z))
#         print("  Not weighed")
# 
#     z_fit = gaussian(_y, A_laser[0], A_laser[1], A_laser[2], A_laser[3])
#     r = _z - z_fit
#     temp2 = numpy.argmax(_z)
# #     print(temp2)
#     temp = (r / weight) ** 2
#     chisq = numpy.sum(temp[(temp2-10):(temp2+10)])
#     
#     print("  Mean:        %5.3f" % A_laser[1])
#     print("  SD:          %5.3f" % A_laser[0])
#     print("  Y-offset:    %5.3f" % A_laser[2])
#     print("  Scale:       %5.3f" % A_laser[3])
#     print("  Chi2:        %5.3f" % chisq)
#     
#     
#     
#     
# 
#     if ax_fit:
#         ax_fit.plot(_y, _z)
#         ax_fit.plot(_y, z_fit)
#     
#     if ax_res:
#         ax_res.plot(_y, r)
# 
#     temp = A_laser[1] - _y[0]
#     x_min = A_laser[1] - 1 
#     if x_min < _y[0]:
#         x_min = _y[0]
#     x_max = A_laser[1] + 1 
#     
#     if ax_fit:
#         ax_fit.set_xlim(x_min, x_max)
#     if ax_res:
#         ax_res.set_xlim(x_min, x_max)
#     
#     if ax_fit and ax_res:
#         if title != "":
#             ax_fit.set_title("Fit of laser of %s" % (title))
#         ax_res.set_xlabel("Time (ns)")
#         ax_fit.set_ylabel("Intensity (counts, averaged)")
#         ax_res.set_ylabel("Residue") 
#     elif ax_fit:
#         if title != "":
#             ax_fit.set_title("Fit of laser of %s" % (title))
#         ax_fit.set_xlabel("Time (ns)")
#         ax_fit.set_ylabel("Intensity (counts, averaged)")
#     elif ax_res:
#         if title != "":
#             ax_res.set_title("Residue of laser of %s" % (title))
#         ax_res.set_xlabel("Time (ns)")
#         ax_res.set_ylabel("Residue") 
# 
#     
#     return A_laser, _y, _z, z_fit




# def gaussian(t, sigma, mu, y, A):
#     """
#     A[0]: sigma (sigma^2 = variance)
#     A[1]: mu (mean)
#     A[2]: offset 
#     A[3]: scale, before offset
# 
#     """
#     res = ( A / (sigma * numpy.sqrt(2*numpy.pi)) ) * numpy.exp( -(t - mu)**2 / (2 * sigma**2) ) + y
#     return res
#     
# def poisson(k, lamb):
#     return (lamb**k/scipy.misc.factorial(k)) * numpy.exp(-lamb)
#     
# 
# def fit_laser_pulse_poisson(x, y, z, xrange = [0,-1], yrange = [0,-1], frame = "data", title = "", ax_fit = False, ax_res = False):
# 
#     _z = get_average(x, y, z, axis = "y", range = xrange, frame = frame)
#     
#     if yrange[0] != 0 and yrange[1] != -1:
#         a, b = CF.find_indices_for_range(y, yrange, frame = frame, round = "maximize")
#         
#         _z = _z[a:b]
#         _y = y[a:b]
#     else:
#         _y = numpy.copy(y)
# 
#     i = numpy.argmax(_z)
#     A_start = [1]
# 
# #     _z += 0.00001
# #     weight = numpy.sqrt(_z)
# #     print(weight == 0)
# #     weight += 0.0000001
#     #     A_start = [1,1,0,1]
#     A_laser, matcov = scipy.optimize.curve_fit(f = poisson, xdata = _y, ydata = _z, p0 = A_start)
#     
#     print(A_laser)
#     
#     z_fit = 0
#     return A_laser, _y, _z, z_fit
# 
# 
# def fit_laser_pulse(x, y, z, xrange = [0,-1], yrange = [0,-1], frame = "data", title = "", ax_fit = False, ax_res = False, flag_weigh_laser = False):
# 
#     _z = get_average(x, y, z, axis = "y", range = xrange, frame = frame)
#     
#     if yrange[0] != 0 and yrange[1] != -1:
#         a, b = CF.find_indices_for_range(y, yrange, frame = frame, round = "maximize")
#         
#         _z = _z[a:b]
#         _y = y[a:b]
#     else:
#         _y = numpy.copy(y)
# 
#     # initial try for laser pulse
#     # SD, mean, offset, scale
#     i = numpy.argmax(_z)
#     A_start = [0.1, _y[i], 0, _z[i]]
# 
#     print("FIT RESULTS LASER %s" % (title))
# 
#     if flag_weigh_laser:
# #         _z += 0.00001
#         
#         weight = numpy.copy(_z)
#         temp = numpy.where(weight == 0)[0]
#         weight[temp] += 1e20
#         weight = numpy.sqrt(weight)
#         A_laser, matcov = scipy.optimize.curve_fit(f = gaussian, xdata = _y, ydata = _z, p0 = A_start, sigma = weight, absolute_sigma = True)
#         print("  Weighed")
#     else:
#         A_laser, matcov = scipy.optimize.curve_fit(f = gaussian, xdata = _y, ydata = _z, p0 = A_start)
#         weight = numpy.ones(len(_z))
#         print("  Not weighed")
# 
#     z_fit = gaussian(_y, A_laser[0], A_laser[1], A_laser[2], A_laser[3])
#     r = _z - z_fit
#     temp2 = numpy.argmax(_z)
# #     print(temp2)
#     temp = (r / weight) ** 2
#     chisq = numpy.sum(temp[(temp2-10):(temp2+10)])
#     
#     print("  Mean:        %5.3f" % A_laser[1])
#     print("  SD:          %5.3f" % A_laser[0])
#     print("  Y-offset:    %5.3f" % A_laser[2])
#     print("  Scale:       %5.3f" % A_laser[3])
#     print("  Chi2:        %5.3f" % chisq)
#     
#     
#     
#     
# 
#     if ax_fit:
#         ax_fit.plot(_y, _z)
#         ax_fit.plot(_y, z_fit)
#     
#     if ax_res:
#         ax_res.plot(_y, r)
# 
#     temp = A_laser[1] - _y[0]
#     x_min = A_laser[1] - 1 
#     if x_min < _y[0]:
#         x_min = _y[0]
#     x_max = A_laser[1] + 1 
#     
#     if ax_fit:
#         ax_fit.set_xlim(x_min, x_max)
#     if ax_res:
#         ax_res.set_xlim(x_min, x_max)
#     
#     if ax_fit and ax_res:
#         if title != "":
#             ax_fit.set_title("Fit of laser of %s" % (title))
#         ax_res.set_xlabel("Time (ns)")
#         ax_fit.set_ylabel("Intensity (counts, averaged)")
#         ax_res.set_ylabel("Residue") 
#     elif ax_fit:
#         if title != "":
#             ax_fit.set_title("Fit of laser of %s" % (title))
#         ax_fit.set_xlabel("Time (ns)")
#         ax_fit.set_ylabel("Intensity (counts, averaged)")
#     elif ax_res:
#         if title != "":
#             ax_res.set_title("Residue of laser of %s" % (title))
#         ax_res.set_xlabel("Time (ns)")
#         ax_res.set_ylabel("Residue") 
# 
#     
#     return A_laser, _y, _z, z_fit
# 
# 
# 
# 
# # def gaussian_single_exponential(t, a, b, c, d, e):
# # 
# #     """
# #     A[0]: sigma (sigma^2 = variance)
# #     A[1]: mu (mean)
# #     A[2]: offset 
# #     A[3]: scale, before offset
# # 
# #     """
# #     
# #     g = numpy.fft.fft(( b / (t[1] * numpy.sqrt(2*numpy.pi)) ) * numpy.exp( -(t[0] - t[2])**2 / (2 * t[1]**2) ) )
# #     e = numpy.fft.fft( d * numpy.exp(-t[0] / a) + c)
# #     
# #     f = g * e
# #     
# #     y = numpy.real(numpy.fft.ifft(f)) 
# # 
# #     return y    
# 
# 
# 
# # def gaussian_single_exponential(t, b, c, d, e, f):
# #     """
# #     t[1] a A[0]: sigma (sigma^2 = variance) 
# #     b A[1]: mu (mean)
# #     c A[2]: offset 
# #     d A[3]: scale, before offset
# #     e: amplitude exponential
# #     f = decay rate
# #     """
# #     g = ( d / (t[1] * numpy.sqrt(2*numpy.pi)) ) * numpy.exp( -(t[0] - b)**2 / (2 * t[1]**2) ) 
# #     h = e * numpy.exp(-t[0] / f) 
# #     _g = numpy.fft.fft(g)
# #     _h = numpy.fft.fft(h)
# #     y = numpy.real(numpy.fft.ifft(_g * _h)) + c
# #     return y
# 
# def gaussian_single_exponential(t, b, c, d, f):
#     """
#     t[1] a A[0]: sigma (sigma^2 = variance) 
#     b A[1]: mu (mean)
#     c A[2]: offset 
#     d A[3]: scale, before offset
#     e: amplitude exponential
#     f = decay rate
#     """
#     g = numpy.exp( -(t[0] - b)**2 / (2 * t[1]**2) ) 
#     h = numpy.exp(-t[0] / f) 
#     _g = numpy.fft.fft(g)
#     _h = numpy.fft.fft(h)
#     y = d * numpy.real(numpy.fft.ifft(_g * _h)) + c
#     return y
# 
# 
# def fit_lifetime(x, y, z, A_laser = [], xrange = [0,-1], yrange = [0,-1], frame = "data"):
# 
#     data = get_average(x, y, z, axis = "y", range = xrange, frame = frame)
#     
#     if yrange[0] != 0 and yrange[1] != -1:
#         a, b = CF.find_indices_for_range(y, yrange, frame = frame, round = "maximize")
#         data = data[a:b]
#         _y = y[a:b]
#     else:
#         _y = y
# 
#     # numpy.savetxt("/Users/rbloem/Transporter/arcnl/Measurements/20170502/time.csv", _y)
# #     numpy.savetxt("/Users/rbloem/Transporter/arcnl/Measurements/20170502/intensity.csv", data)
# 
#     t = [_y, A_laser[0]] #, A_laser[2]]
#     A_start = scipy.array([A_laser[1], A_laser[2], A_laser[3], y[-1]/5, ])
#     print(A_start)
#     
#     sigma = 1/numpy.sqrt(data)
# 
# #     data -= numpy.amin(data)
# #     data /= numpy.amax(data)
# 
# #     temp = numpy.isnan(sigma)
# #     print(temp)
# 
#     A_out, matcov = scipy.optimize.curve_fit(gaussian_single_exponential, t, data, p0 = A_start, sigma = sigma, absolute_sigma = False)
# 
# #     print("sigma       {:4.2f}".format(A_out[0]))
# #     print("mean        {:4.2f}".format(A_out[1]))
# #     print("offset      {:4.2f}".format(A_out[2]))
# #     print("scale       {:4.2f}".format(A_out[3]))
# #     print("amplitude   {:4.2f}".format(A_out[4]))
#     print("decay rate  {:4.2f}".format(1/A_out[3]))
#     print("lifetime    {:4.2f}".format(A_out[3]))
#     
#  #        t[1] a A[0]: sigma (sigma^2 = variance) 
# #     b A[1]: mu (mean)
# #     c A[2]: offset 
# #     d A[3]: scale, before offset
# #     e: amplitude exponential
# #     f = decay rate
# 
#     
#     
#     print("single exp gaus", A_out)
#     out = gaussian_single_exponential(t, A_out[0], A_out[1], A_out[2], A_out[3])
# 
#     r = data - out
# 
#     fig = plt.figure()
#     
#     ax = [0] * 2
#     ax[0] = fig.add_subplot(211)
#     ax[1] = fig.add_subplot(212)
#     
#     ax[0].plot(_y, data)
#     ax[0].plot(_y, out)
#     
#     xlim = ax[0].get_xlim()
#     ylim = ax[0].get_ylim()
#     
#     x = (xlim[0] + xlim[1]) /2
#     y = (ylim[0] + ylim[1]) /2
#     
#     s = "lifetime {:4.2f} ns".format(A_out[3])
#     ax[0].text(x,y, s = s)
#     
#     ax[1].plot(_y, r)
#     
#     return fig, ax, _y, out
# 
# 
# #     plt.plot(_y, data)
# #     plt.plot(_y, out)


if __name__ == "__main__": 
    pass