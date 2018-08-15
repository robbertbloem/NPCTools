"""

"""

import struct
import re

import scipy
import numpy
import matplotlib 
import matplotlib.pyplot as plt

import NPCTools.Resources.CommonFunctions as CF


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
    

        


def import_streak(paf):
    """
    Import the raw data and extract all information. 
    """

    # read file
    with open(paf, 'rb') as f:
        raw = f.read()
    
    comment_length = struct.unpack_from("h", raw, offset = 2)[0] 
    width = struct.unpack_from("h", raw, offset = 4)[0] 
    height = struct.unpack_from("h", raw, offset = 6)[0] 
    wl_offset = struct.unpack_from("h", raw, offset = 8)[0] 
    t_offset = struct.unpack_from("h", raw, offset = 10)[0] 
    type = struct.unpack_from("h", raw, offset = 12)[0] 

    comment_string = str(raw[64:comment_length+64])
    meta_dict = parse_comment(comment_string)

    meta_dict["Processing"]["Filename"] = paf

    data = struct.unpack_from("i" * width * height, raw, offset = comment_length+64)
    data = numpy.asarray(data, dtype = numpy.float64) # save as float
    data = numpy.reshape(data, (height, width)) # reshape from 1D to 2D

    w_axis_offset = (comment_length + 64) + width * height * 4
    t_axis_offset = w_axis_offset + width * 4

    w_axis = struct.unpack_from("f" * width, raw, offset = w_axis_offset) 
    t_axis = struct.unpack_from("f" * height, raw, offset = t_axis_offset)

    w_axis = numpy.array(w_axis)
    t_axis = numpy.array(t_axis)
    
    if w_axis[0] > w_axis[-1]:
        w_axis = w_axis[::-1]
        data = data[:,::-1]

    return data, w_axis, t_axis, meta_dict



def get_area(w_axis, t_axis, z, w_range = [0,-1], t_range = [0,-1], frame = "data", debug = 0):
    """
    w_axis:
    t_axis, z
    w_range: [min, max]. [0,-1] means the complete range.
    t_range = [min, max]. [0,-1] means the complete range.
    frame = "z"
    """

    if debug > 0:
        print("NPCTools.Resources.StreakFunctions.get_area") 

    if w_range[0] == 0 and w_range[1] == -1:
        _z = z[:,:]
        _w_axis = w_axis[:]
    else:
        a, b = CF.find_indices_for_range(w_axis, w_range, frame = frame, round = "maximize")    
        if debug > 0:
            print("  w_axis {:.1f} nm [{:d}] to {:.1f} nm [{:d}]".format( w_axis[a], a, w_axis[b], b)) 
        _z = z[:, a:b]
        _w_axis = w_axis[a:b]

    if t_range[0] == 0 and t_range[1] == -1:
        _z = _z[:,:]
        _t_axis = t_axis[:]
    else:
        a, b = CF.find_indices_for_range(t_axis, t_range, frame = frame, round = "maximize")
        if debug > 0:
            print("  t_axis time from {:.1f} [{:d}] to {:.1f} [{:d}] ".format(t_axis[a], a, t_axis[b], b)) 
        _z = _z[a:b, :]
        _t_axis = t_axis[a:b]

    return _w_axis, _t_axis, _z




def get_average(w_axis, t_axis, z, axis = "w", range = [0,-1], frame = "data", return_cumulative = False, debug = 0):
    """
    INPUTS:
    w_axis: wavelength axis
    t_axis: time axis
    z: data
    axis: 'w' or 't', the axis you want to keep. Average w means that the t-axis is averaged. 
    range: start and end of the range to be averaged.
    frame: 'data' or 'index'
    """

    if debug > 0:
        print("NPCTools.Resources.StreakFunctions.get_average") 

    if axis in ["w", "w_axis"]:
        if range[0] == 0 and range[1] == -1:
            _z = numpy.mean(z, axis = 0)
            _s = numpy.sum(z, axis = 0)
        else:
            a, b = CF.find_indices_for_range(t_axis, range, frame = frame, round = "maximize")   
            if debug > 0:
                print("  w_axis: from wavelength {:4.1f} [{:d}] to {:4.1f} [{:d}] ".format(w_axis[a], a, w_axis[b], b)) 
            _z = numpy.mean(z[a:b, :], axis = 0)
            _s = numpy.sum(z[a:b, :], axis = 0)
    else:
        if range[0] == 0 and range[1] == -1:
            _z = numpy.mean(z, axis = 1)
            _s = numpy.sum(z, axis = 1)
        else:
            a, b = CF.find_indices_for_range(w_axis, range, frame = frame, round = "maximize")
            if debug > 0:
                print("  t_axis time from {:.1f} [{:d}] to {:.1f} [{:d}] ".format(t_axis[a], a, t_axis[b], b)) 
            _z = numpy.mean(z[:, a:b], axis = 1)
            _s = numpy.sum(z[:, a:b], axis = 1)
    if return_cumulative:
        return _s
    else:
        return _z



def gaussian_single_exponential(t, a, b, c, d, e, f):
    """
    t[0]: time axis
    Laser pulse
    a A[0]: sigma (sigma^2 = variance) 
    b A[1]: mu (mean)
    c A[2]: offset 
    d A[3]: scale, before offset
    
    e: amplitude exp 1
    f: decay rate exp 1
    """
    
    return gaussian_single_exponential_helper(t, a, b, d, e, f) + c


def gaussian_single_exponential_helper(t, a, b, d, e, f):
    """
    t[0]: time axis
    Laser pulse
    a A[0]: sigma (sigma^2 = variance) 
    b A[1]: mu (mean)
    c A[2]: offset 
    d A[3]: scale, before offset
    
    e: amplitude exp 1
    f: decay rate exp 1
    """
    g = numpy.exp( -(t[0] - b)**2 / (2 * a**2) ) 
    h = e * numpy.exp(-t[0] / f) 
    _g = numpy.fft.fft(g)
    _h = numpy.fft.fft(h)
    y = d * numpy.real(numpy.fft.ifft(_g * _h))
    return y
    


def gaussian_double_exponential(t, a, b, c, d, e, f, k, m):
    """
    t[0]: time axis
    Laser pulse
    a A[0]: sigma (sigma^2 = variance) 
    b A[1]: mu (mean)
    c A[2]: offset 
    d A[3]: scale, before offset
    
    e: amplitude exp 1
    f: decay rate exp 1
    
    k: amplitude exp 2
    m: decay rate exp 2
    """
    
    return gaussian_single_exponential_helper(t, a, b, d, e, f) + gaussian_single_exponential_helper(t, a, b, d, k, m) + c


def gaussian_triple_exponential(t, a, b, c, d, e, f, k, m, n, p):
    """
    t[0]: time axis
    Laser pulse
    a A[0]: sigma (sigma^2 = variance) 
    b A[1]: mu (mean)
    c A[2]: offset 
    d A[3]: scale, before offset
    
    e: amplitude exp 1
    f: decay rate exp 1
    
    k: amplitude exp 2
    m: decay rate exp 2

    n: amplitude exp 3
    p: decay rate exp 3
    """
    
    return gaussian_single_exponential_helper(t, a, b, d, e, f) + gaussian_single_exponential_helper(t, a, b, d, k, m) + gaussian_single_exponential_helper(t, a, b, d, n, p) + c






def fit_lifetime(w_axis, t_axis, z, w_range = [0,-1], t_range = [0,-1], A_laser = [], frame = "data", flag_weigh_laser = False, debug = 0, fit_type = "single_exp", exp_hints = -1):
    
    if debug > 0:
        print("NPCTools.Resources.StreakFunctions.fit_lifetime") 
        print("  fit type: {:s}".format(fit_type))

    _w_axis, _t_axis, _z = get_area(w_axis, t_axis, z, w_range = w_range, t_range = t_range, frame = frame, debug = debug)

    _z = get_average(_w_axis, _t_axis, _z, axis = "t", range = [0,-1], frame = frame, debug = debug)
    
    if len(A_laser) == 0:       
        A_laser = [0] * 4
        A_laser[0] = 0.1  # sigma
        A_laser[1] = _t_axis[numpy.argmax(_z)] # mean
        A_laser[2] = 0 # offset
        A_laser[3] = 1 # scale
        

    
    sigma = 1/numpy.sqrt(_z)
    
    if fit_type == "single_exp":
        t = [_t_axis]
        A_start = scipy.array([A_laser[0], A_laser[1], A_laser[2], A_laser[3], numpy.amax(_z), _t_axis[-1]/5, ])
        A_out, matcov = scipy.optimize.curve_fit(gaussian_single_exponential, t, _z, p0 = A_start, sigma = sigma, absolute_sigma = False)
        
        y_fit = gaussian_single_exponential(t, A_out[0], A_out[1], A_out[2], A_out[3], A_out[4], A_out[5])
        
        print("Laser")
        print("  Sigma:     {:3.2f}".format(A_out[0]))
        print("  Mean:      {:3.2f}".format(A_out[1]))
        print("  Offset:    {:3.2f}".format(A_out[2]))
        print("  Scale:     {:3.2f}".format(A_out[3]))
        print("Exponential")
        print("  Amplitude: {:3.2f}".format(A_out[4]))
        print("  Lifetime:  {:3.2f}".format(A_out[5]))
        
    
    
    elif fit_type == "double_exp":
        t = [_t_axis]
        
        if type(exp_hints) == list:
            hint1 = exp_hints[0]
            hint2 = exp_hints[1]
        else:
            hint1 = _t_axis[-1]/5
            hint2 = _t_axis[-1]/20
        
        A_start = scipy.array([A_laser[0], A_laser[1], A_laser[2], A_laser[3], numpy.amax(_z), hint1, numpy.amax(_z)/2, hint2,])
        A_out, matcov = scipy.optimize.curve_fit(gaussian_double_exponential, t, _z, p0 = A_start, sigma = sigma, absolute_sigma = False)
        
        y_fit = gaussian_double_exponential(t, A_out[0], A_out[1], A_out[2], A_out[3], A_out[4], A_out[5], A_out[6], A_out[7])


        print("Laser")
        print("  Sigma:     {:3.2f}".format(A_out[0]))
        print("  Mean:      {:3.2f}".format(A_out[1]))
        print("  Offset:    {:3.2f}".format(A_out[2]))
        print("  Amplitude: {:3.2f}".format(A_out[3]))
        print("Exponential 1")
        print("  Amplitude: {:3.2f}".format(A_out[4]))
        print("  Lifetime:  {:3.2f}".format(A_out[5]))
        print("Exponential 2")
        print("  Amplitude: {:3.2f}".format(A_out[6]))
        print("  Lifetime:  {:3.2f}".format(A_out[7]))


    elif fit_type == "triple_exp":
        t = [_t_axis]
        
        if type(exp_hints) == list:
            hint1 = exp_hints[0]
            hint2 = exp_hints[1]
            hint2 = exp_hints[2]
        else:
            hint1 = _t_axis[-1]/5
            hint2 = _t_axis[-1]/10
            hint2 = _t_axis[-1]/20
        
        A_start = scipy.array([A_laser[0], A_laser[1], A_laser[2], A_laser[3], numpy.amax(_z), hint1, numpy.amax(_z)/2, hint2, numpy.amax(_z)/2, hint3,])
        
        A_out, matcov = scipy.optimize.curve_fit(gaussian_double_exponential, t, _z, p0 = A_start, sigma = sigma, absolute_sigma = False)
        
        y_fit = gaussian_double_exponential(t, A_out[0], A_out[1], A_out[2], A_out[3], A_out[4], A_out[5], A_out[6], A_out[7], A_out[8], A_out[9])


        print("Laser")
        print("  Sigma:     {:3.2f}".format(A_out[0]))
        print("  Mean:      {:3.2f}".format(A_out[1]))
        print("  Offset:    {:3.2f}".format(A_out[2]))
        print("  Amplitude: {:3.2f}".format(A_out[3]))
        print("Exponential 1")
        print("  Amplitude: {:3.2f}".format(A_out[4]))
        print("  Lifetime:  {:3.2f}".format(A_out[5]))
        print("Exponential 2")
        print("  Amplitude: {:3.2f}".format(A_out[6]))
        print("  Lifetime:  {:3.2f}".format(A_out[7]))
        print("Exponential 3")
        print("  Amplitude: {:3.2f}".format(A_out[8]))
        print("  Lifetime:  {:3.2f}".format(A_out[9]))



    
    return _t_axis, _z, y_fit, A_out





if __name__ == "__main__": 
    pass