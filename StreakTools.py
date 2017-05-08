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


imp.reload(CF)

def parse_comment(string):
    """
    The comment in the streak files is a long string. In this function the string is separated into sections with a header. For each header the keys and values are found using regex. 
    The output meta_dict is a dictionary with the headers as keys. The value is another dict which contains the key and values. 
    
    """
    
    p = re.compile(",[a-zA-Z\s\.]*=")
    
    headers = ["Application", "Camera", "Acquisition", "Grabber", "DisplayLUT", "ExternalDevices", "Streak camera", "Spectrograph", "Delay box", "Delay2 box", "Scaling", "Comment"] 
    headers = headers[::-1]
    
    meta_dict = {}

    for h in range(len(headers)):
        splitstring = string.split("[%s]" % headers[h])
        string = splitstring[0]
        headerstring = splitstring[1]
        
        temp_dict = {}

        regexresult = p.findall(headerstring)
        for r in range(len(regexresult)):
            regexresult[r] = regexresult[r][1:-1]
    
        for r in range(len(regexresult)):
            x = re.findall(',%s=\"[a-zA-Z\s\.,:\-0-9/#\(\)\_]*\"' % regexresult[r], headerstring)
            if len(x) == 0:
                x = re.findall(',%s=[a-zA-Z0-9#\-\(\)\_]*' % regexresult[r], headerstring)
                
            a, b = x[0].split("=")
            if len(b) > 1:
                if b[0] == '"':                
                    b = b[1:-1]
            
            temp_dict[regexresult[r]] = b
            
        meta_dict[headers[h]] = temp_dict
        
    meta_dict["Processing"] = {}

    return meta_dict
    
def print_meta_dict(meta_dict):
    
    for k, v in meta_dict.items():
        print("%s" % k.upper())
        for x, y in v.items():
            print("\t%25s = %s" % (x, y))
        




def import_streak(paf):

    with open(paf, 'rb') as f:
        data = f.read()

    comment_length = struct.unpack_from("h", data, offset = 2)[0]
    width = struct.unpack_from("h", data, offset = 4)[0]
    height = struct.unpack_from("h", data, offset = 6)[0]
    x_offset = struct.unpack_from("h", data, offset = 8)[0]
    y_offset = struct.unpack_from("h", data, offset = 10)[0]
    type = struct.unpack_from("h", data, offset = 12)[0]

#     print(comment_length, width, height, x_offset, y_offset, type)

    comment_string = str(data[64:comment_length+64])
    meta_dict = parse_comment(comment_string)

    meta_dict["Processing"]["Filename"] = paf

    data2 = struct.unpack_from("i" * width * height, data, offset = comment_length+64)

    data2 = numpy.reshape(data2, (height, width))

    x_axis_offset = (comment_length + 64) + width * height * 4
    y_axis_offset = x_axis_offset + width * 4

    x_axis = struct.unpack_from("f" * width, data, offset = x_axis_offset) 
    y_axis = struct.unpack_from("f" * height, data, offset = y_axis_offset)

    x_axis = numpy.array(x_axis)
    y_axis = numpy.array(y_axis)
    
    if x_axis[0] > x_axis[-1]:
        x_axis = x_axis[::-1]
        data2 = data2[:,::-1]

    return data2, x_axis, y_axis, meta_dict






def get_average(x, y, z, axis = "x", range = [0,-1], frame = "data"):
    """
    INPUTS:
    x: x-axis
    y: y-axis
    z: data
    axis: 'x' or 'y', the axis you want to keep. Average x means that the y-axis is averaged. 
    range: start and end of the range to be averaged.
    frame: 'data' or 'index'
    """


    if axis == "x":
        a, b = CF.find_indices_for_range(y, range, frame = frame, round = "maximize")    
        data = numpy.mean(z[a:b, :], axis = 0)
    else:
        a, b = CF.find_indices_for_range(x, range, frame = frame, round = "maximize")
        data = numpy.mean(z[:, a:b], axis = 1)
        
    return data
    



def plot_streak(x, y, z, meta_dict, ax = False, **kwargs):

    """
    x_range, y_range (list with 2 elements, default [None, None]): 
    
    """

    if ax == False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    
    if "z_range" in kwargs:
        levels = numpy.linspace(kwargs["z_range"][0], kwargs["z_range"][1])
        kwargs["levels"] = levels
        
    ax.contourf(x, y, z, **kwargs)
    
    ax.set_xlabel("Wavelength (%s)" % meta_dict["Scaling"]["ScalingXUnit"])
    ax.set_ylabel("Time (%s)" % meta_dict["Scaling"]["ScalingYUnit"])

    if "x_range" in kwargs:
        ax.set_xlim(kwargs["x_range"])

    if "y_range" in kwargs:
        ax.set_ylim(kwargs["y_range"])

    return ax


def fit_laser_pulse(x, y, z, xrange = [0,-1], yrange = [0,-1], frame = "data"):

    data = get_average(x, y, z, axis = "y", range = xrange, frame = frame)
    
    if yrange[0] != 0 and yrange[1] != -1:
        a, b = CF.find_indices_for_range(y, yrange, frame = frame, round = "maximize")
        
        data = data[a:b]
        _y = y[a:b]
    else:
        _y = y
    
    # initial try for laser pulse
    # SD, mean, offset, scale
    i = numpy.argmax(data)
    A_start = [0.1, _y[i], 0, 1]

    # fit
    A_laser = MATH.fit(_y, data, EQ.rb_gaussian, A_start, return_all = False)
    
    # check output
    out = EQ.rb_gaussian(A_laser, _y)
    
    print("FIT RESULTS")
    print("Mean:", A_laser[1])
    print("SD:", A_laser[0])
    print("Offset:", A_laser[2])
    print("Scale:", A_laser[3])
    
    plt.plot(_y, data)
    plt.plot(_y, out)
    
    return A_laser




def gaussian_single_exponential(t, a, b, c):

    """
    A[0]: sigma (sigma^2 = variance)
    A[1]: mu (mean)
    A[2]: offset 
    A[3]: scale, before offset

    """
    
    g = numpy.fft.fft(( b / (t[1] * numpy.sqrt(2*numpy.pi)) ) * numpy.exp( -(t[0] - t[2])**2 / (2 * t[1]**2) ) )
    e = numpy.fft.fft( numpy.exp(-t[0] / a))
    
    f = g * e
    
    y = numpy.real(numpy.fft.ifft(f)) + c

    return y    


def fit_lifetime(x, y, z, A_laser = [], xrange = [0,-1], yrange = [0,-1], frame = "data"):

    data = get_average(x, y, z, axis = "y", range = xrange, frame = frame)
    
    if yrange[0] != 0 and yrange[1] != -1:
        a, b = CF.find_indices_for_range(y, yrange, frame = frame, round = "maximize")
        data = data[a:b]
        _y = y[a:b]
    else:
        _y = y

    t = [_y, A_laser[0], A_laser[1]] #, A_laser[2]]
    A = scipy.array([y[-1]/4, A_laser[3], 0])
    print(A)
    param = (t, data, gaussian_single_exponential)
    
#     A_out, cov_x, infodict, mesg, ier = scipy.optimize.leastsq(MATH.minimize, A, args=param, full_output=True)
    
    
    A_out, matcov = scipy.optimize.curve_fit(gaussian_single_exponential, t, data, p0 = A)
    
    
    
#     A_out = scipy.optimize.leastsq(scipy.optimize.minimize, A, args = param)
    
#     A_out = scipy.optimize.curve_fit(gaussian_single_exponential, t, data, p0 = A)[0]



    print(len(A_out))
    print(A_out)
    out = gaussian_single_exponential(t, A_out[0], A_out[1], A_out[2])

    plt.plot(_y, data)
    plt.plot(_y, out)


if __name__ == "__main__": 
    pass