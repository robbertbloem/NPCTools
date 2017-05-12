from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

import NPCTools.Resources.Mathematics as MATH
import NPCTools.Resources.Equations as EQ

def gaussian_beam_diameter(position, intensity, ax = False, label = "", flag_interpolate = True, interpolate_kind = "linear", interpolate_step = -1, position_unit = "mm", center_plot_on_zero = True): 
    """
    Calculate the diameter of a gaussian beam.
    
    INPUT:
    - position (list, ndarray): list with positions
    - intensity (list, ndarray): list with intensities corresponding to the positions
    - ax (axis object, or False): if False, the result will not be plotted
    - label (string): label used in the fit
    - flag_interpolate (BOOL, True): if the positions are not equidistant, then it has to be interpolated. 
    - interpolate_kind
    - interpolate_step (number): stepsize of the interpolation. If the value is negative, that is the number of steps that will be used. If the value is 0, 100 steps will be used. 
    - position_unit (string, mm): unit of the position
    - center_plot_on_zero (BOOL, True): Subtract the mean position from the position. 
    
    
    """
    
    print("--- %s ---" % label)
    
    if position[1] < position[0]:
        position = position[::-1]
        intensity = intensity[::-1]
    

    if flag_interpolate:
        if interpolate_step == 0:
            x = numpy.linspace(position[0], position[-1], num = 100)
        elif interpolate_step < 0:
            x = numpy.linspace(position[0], position[-1], num = -interpolate_step)
        else:
            x = numpy.arange(position[0], position[-1], interpolate_step)
        f = interp1d(position, intensity, kind = interpolate_kind)
        y = f(x)
    else:
        x = numpy.copy(position)
        y = numpy.copy(intensity)
        
    dx = x[1] - x[0]

    mean_guess = x[numpy.argmax(y)]
    A = [1, mean_guess, 1, 1]
    A_final = MATH.fit(x, y, EQ.rb_gaussian, A)

    
    
    if center_plot_on_zero:
        print("Mean: %2.3f %s (shifted in plot)" % (A_final[1], position_unit))
    else:
        print("Mean: %2.3f %s" % (A_final[1], position_unit))
    
    y_fit = EQ.rb_gaussian(A_final, x)
    
    if center_plot_on_zero:
        x -= A_final[1]
        A_final[1] = 0
    
    temp = numpy.where(y_fit > (0.135 * numpy.amax(y_fit)))[0]
    print("1/e2: %2.3f %s" % ( x[temp[-1]] - x[temp[0]], position_unit))
    
    temp = numpy.where(y_fit > (0.5 * numpy.amax(y_fit)))[0]
    print("FWHM: %2.3f %s" % (x[temp[-1]] - x[temp[0]], position_unit))    

    y_fit -= numpy.amin(y_fit)
    y /= numpy.amax(y_fit)
    y_fit /= numpy.amax(y_fit)
    
    s = label
    ax.plot(x, y, label = s)
    s = label + " fit"
    ax.plot(x, y_fit, label = s)
    
    ax.set_xlabel("Position (%s, shifted)" % position_unit)
    ax.set_ylabel("Intensity (norm)")
    
    ax.legend()
    
    print("Stepsize positions: %2.3f" % dx) 
    
    print("Fitting parameters:")
    print("    Sigma:", A_final[0])
    print("    Mean:", A_final[1])
    print("    Y-offset:", A_final[2])
    print("    Scale:", A_final[3])
    

#     
#     print("mean:", (x_offset + A_final[1])) 
#     
#     print(A_final)  


if __name__ == "__main__": 
    pass