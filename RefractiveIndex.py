from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

import NPCTools.Resources.RI_read_yaml as RIRY
import NPCTools.Resources.Constants as CONST
import NPCTools.Resources.CommonFunctions as CF
import NPCTools.Debug as DEBUG

def interpolate_data(db_record, wl_um, interpolate_kind = "default", verbose = 0):
    """
    db_record contains the original data
    wl_um are the wavelengths we want to know
    kind = linear, quadratic, cubic
    """    
    
    if interpolate_kind == "default":
        interpolate_kind = "linear"
    
    DEBUG.verbose("  Interpolating data using %s" % (interpolate_kind), verbose_level = 1)
    
    f = interp1d(db_record["data"][:,0], db_record["data"][:,1], kind = interpolate_kind)
    ri = f(wl_um)    

    return ri



def ri_for_wavelengths(db_record, wl_um, interpolate_kind = "default", verbose = 0):
    """
    db_record is the dictionary
    """

    if "type" not in db_record:
        raise KeyError("ri_for_wavelengths(): The database record does not have a key 'type'.")
        
    wl_um = CF.make_numpy_ndarray(wl_um)
   
    error_string = "Calculation of the refractive index is not implemented for "
    
    # check the range 
    if "formula" in db_record["type"]:
        
        if wl_um[0] < db_record["range"][0]: 
            raise ValueError("Error, wavelength %1.2f micron is too low! It should be above %1.2f micron." % (wl_um[0], db_record["range"][0]) )
        elif wl_um[-1] > db_record["range"][1]:
            raise ValueError("Error, wavelength %1.2f micron is too high! It should be below %1.2f micron." % (wl_um[-1], db_record["range"][-1]) )
  
    if "tabulated" in db_record["type"]:
        if wl_um[0] < db_record["data"][:,0][0]: 
            raise ValueError("Error, wavelength %1.2f micron is too low! It should be above %1.2f micron." % (wl_um[0], db_record["data"][:,0][0]) )
        elif wl_um[-1] > db_record["data"][:,0][-1]:
            raise ValueError("Error, wavelength %1.2f micron is too high! It should be below %1.2f micron." % (wl_um[-1], db_record["data"][:,0][-1]) )   

    if db_record["type"] == "formula 1":
        DEBUG.verbose("  Using formula 1 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.ones(len(wl_um)) + db_record["coefficients"][0] 
        l = wl_um**2
        for i in range(n_terms):
            ri += (db_record["coefficients"][2*i+1] * l) / (l - db_record["coefficients"][2*i+2]**2)
        ri = numpy.sqrt(ri)
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])

    elif db_record["type"] == "formula 2":
        DEBUG.verbose("  Using formula 2 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.ones(len(wl_um)) + db_record["coefficients"][0] 
        l = wl_um**2
        for i in range(n_terms):
            ri += (db_record["coefficients"][2*i+1] * l) / (l - db_record["coefficients"][2*i+2])
        ri = numpy.sqrt(ri)

        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])

    elif db_record["type"] == "formula 3":
        DEBUG.verbose("  Using formula 3 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.zeros(len(wl_um)) + db_record["coefficients"][0] 
        for i in range(n_terms):
            ri += db_record["coefficients"][2*i+1] * wl_um**db_record["coefficients"][2*i+2]
        ri = numpy.sqrt(ri)
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])
    
    elif db_record["type"] == "formula 4":  
        DEBUG.verbose("  Using formula 4 to calculate refractive indices", verbose_level = 1)
        ri = numpy.zeros(len(wl_um)) + db_record["coefficients"][0] 
        n_coeff = len(db_record["coefficients"])
        if n_coeff >= 5:
            ri += (db_record["coefficients"][1] * wl_um**db_record["coefficients"][2]) / ( wl_um**2  - db_record["coefficients"][3]**db_record["coefficients"][4])
        if n_coeff >= 9:
            ri += (db_record["coefficients"][5] * wl_um**db_record["coefficients"][6]) / ( wl_um**2  - db_record["coefficients"][7]**db_record["coefficients"][8])
        n_terms = int((len(db_record["coefficients"]) - 9)/ 2)        
        for i in range(n_terms):
            ri += db_record["coefficients"][2*i+9] * wl_um**db_record["coefficients"][2*i+10]
        ri = numpy.sqrt(ri)
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])
    
    elif db_record["type"] == "formula 5":
        DEBUG.verbose("  Using formula 5 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.zeros(len(wl_um)) + db_record["coefficients"][0] 
        for i in range(n_terms):
            ri += db_record["coefficients"][2*i+1] * wl_um**db_record["coefficients"][2*i+2]
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])
        
    elif db_record["type"] == "formula 6":
        DEBUG.verbose("  Using formula 6 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.ones(len(wl_um)) + db_record["coefficients"][0] 
        l = wl_um**-2
        for i in range(n_terms):
            ri += (db_record["coefficients"][2*i+1]) / (db_record["coefficients"][2*i+2] - l)
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])
        
    elif db_record["type"] == "formula 7":
        raise NotImplementedError(error_string + "formula 7")  
        
    elif db_record["type"] == "formula 8":
        raise NotImplementedError(error_string + "formula 8")  
        
    elif db_record["type"] == "formula 9":
        raise NotImplementedError(error_string + "formula 9")  

    elif db_record["type"] == "tabulated n":
        DEBUG.verbose("  Tabulated data (type n)", verbose_level = 1)
        ri = interpolate_data(db_record, wl_um)
        
    elif db_record["type"] == "tabulated nk":
        DEBUG.verbose("  Tabulated data (type nk)", verbose_level = 1)
        ri = interpolate_data(db_record, wl_um)
        
    else:
        raise ValueError("The type of data in the database record is unknown (usually formula 1-9 or tabulated data). Type here is %s." % db_record["type"]) 
    
    return ri


def gvd_formula_1(x, s): 
    """
    Formula 1 and 2 are the same, except that some coefficients are squared. 
    """
    if len(s) < 9:
        _s = numpy.zeros(9)
        _s[:len(s)] = s
        s = _s

    y = (
        (6 * s[1] * s[2]**2) / (x**4 * (1 - s[2]**2 / x**2)**2) 
        + (8 * s[1] * s[2]**4) / (x**6 * (1 - s[2]**2 / x**2)**3) 
        + (6 * s[3] * s[4]**2) / (x**4 * (1 - s[4]**2 / x**2)**2) 
        + (8 * s[3] * s[4]**4) / (x**6 * (1 - s[4]**2 / x**2)**3) 
        + (6 * s[5] * s[6]**2) / (x**4 * (1 - s[6]**2 / x**2)**2) 
        + (8 * s[5] * s[6]**4) / (x**6 * (1 - s[6]**2 / x**2)**3) 
        + (6 * s[7] * s[8]**2) / (x**4 * (1 - s[8]**2 / x**2)**2) 
        + (8 * s[7] * s[8]**4) / (x**6 * (1 - s[8]**2 / x**2)**3)
    ) / (
        2 * numpy.sqrt(
            s[1] / (1 - s[2]**2 / x**2) 
            + s[3] / (1 - s[4]**2 / x**2) 
            + s[5] / (1 - s[6]**2 / x**2) 
            + s[7] / (1 - s[8]**2 / x**2) 
            + s[0] 
            + 1
        )
    ) - (
        -(2 * s[1] * s[2]**2)/(x**3 * (1 - s[2]**2 / x**2)**2) 
        - (2 * s[3] * s[4]**2)/(x**3 * (1 - s[4]**2 / x**2)**2) 
        - (2 * s[5] * s[6]**2)/(x**3 * (1 - s[6]**2 / x**2)**2) 
        - (2 * s[7] * s[8]**2)/(x**3 * (1 - s[8]**2 / x**2)**2)
    )**2 / (
        4 * (
            s[1]/(1 - s[2]**2 / x**2) 
            + s[3] / (1 - s[4]**2 / x**2) 
            + s[5] / (1 - s[6]**2 / x**2) 
            + s[7] / (1 - s[8]**2 / x**2) 
            + s[0] 
            + 1
        )**(3/2)
    )
    
    return y


def gvd_formula_3(x, s):

    """
    http://www.wolframalpha.com/input/?i=second+derivative+of+sqrt(a+%2B+b*x%5Ec+%2B+d*x%5Ee+%2B+f*x%5Eg+%2B+h*x%5Ei+%2B+j*x%5Ek+%2B+l*x%5Em+%2B+n*x%5Eo+%2Bp*x%5Eq)+with+respect+to+x
    
    """

    if len(s) < 17:
        _s = numpy.zeros(17, dtype = "float")
        _s[:len(s)] = s
        s = _s

    y = (
            s[1] * (s[2] - 1) * s[2] * x**(s[2] - 2) 
            + s[3] * (s[4] - 1) * s[4] * x**(s[4] - 2) 
            + s[5] * (s[6] - 1) * s[6] * x**(s[6] - 2) 
            + s[7] * (s[8] - 1) * s[8] * x**(s[8] - 2) 
            + s[9] * (s[10] - 1) * s[10] * x**(s[10] - 2) 
            + s[11] * (s[12] - 1) * s[12] * x**(s[12] - 2) 
            + s[13] * (s[14] - 1) * s[14] * x**(s[14] - 2) 
            + s[15] * (s[16] - 1) * s[16] * x**(s[16] - 2)
        ) / (
            2 * numpy.sqrt(
                s[0] 
                + s[1] * x**s[2] 
                + s[3] * x**s[4] 
                + s[5] * x**s[6] 
                + s[7] * x**s[8] 
                + s[9] * x**s[10] 
                + s[11] * x**s[12] 
                + s[13] * x**s[14] 
                + s[15] * x**s[16]
            )
        ) - (
            s[1] * s[2] * x**(s[2] - 1) 
            + s[3] * s[4] * x**(s[4] - 1) 
            + s[5] * s[6] * x**(s[6] - 1) 
            + s[7] * s[8] * x**(s[8] - 1) 
            + s[9] * s[10] * x**(s[10] - 1) 
            + s[11] * s[12] * x**(s[12] - 1) 
            + s[13] * s[14] * x**(s[14] - 1) 
            + s[15] * s[16] * x**(s[16] - 1)
        )**2 / (
            4 * (
                s[0] 
                + s[1] * x**s[2] 
                + s[3] * x**s[4] 
                + s[5] * x**s[6] 
                + s[7] * x**s[8] 
                + s[9] * x**s[10] 
                + s[11] * x**s[12] 
                + s[13] * x**s[14] 
                + s[15] * x**s[16]
            )**(3/2)
        )
    return y


def gvd_formula_4(x, s):
    """
    http://www.wolframalpha.com/input/?i=second+derivative+of+sqrt(a%2B+(b*x%5Ec)%2F(x%5E2-d%5Ee)+%2B+(f*x%5Eg)%2F(x%5E2-h%5Ei)+%2B+j*x%5Ek+%2B+l*x%5Em+%2B+n*x%5Eo+%2B+p*x%5Eq+)+with+respect+to+x
    
    
    """
    if len(s) < 17:
        _s = numpy.zeros(17, dtype = "float")
        _s[:len(s)] = s
        s = _s

    y = (
            (s[1] * (s[2] - 1) * s[2] * x**(s[2] - 2))/(x**2 - s[3]**s[4]) 
            + (8 * s[1] * x**(s[2] + 2))/(x**2 - s[3]**s[4])**3 
            - (2 * s[1] * s[2] * x**s[2])/(x**2 - s[3]**s[4])**2 
            - (2 * s[1] * (s[2] + 1) * x**s[2])/(x**2 - s[3]**s[4])**2 
            + (s[5] * (s[6] - 1) * s[6] * x**(s[6] - 2))/(x**2 - s[7]**s[8]) 
            + (8 * s[5] * x**(s[6] + 2))/(x**2 - s[7]**s[8])**3 
            - (2 * s[5] * s[6] * x**s[6])/(x**2 - s[7]**s[8])**2 
            - (2 * s[5] * (s[6] + 1) * x**s[6])/(x**2 - s[7]**s[8])**2 
            + s[9] * (s[10] - 1) * s[10] * x**(s[10] - 2) 
            + s[11] * (s[12] - 1) * s[12] * x**(s[12] - 2) 
            + s[13] * (s[14] - 1) * s[14] * x**(s[14] - 2) 
            + s[15] * (s[16] - 1) * s[16] * x**(s[16] - 2)
        ) / (
            2 * numpy.sqrt(
                s[0] 
                + (s[1] * x**s[2])/(x**2 - s[3]**s[4]) 
                + (s[5] * x**s[6])/(x**2 - s[7]**s[8]) 
                + s[9] * x**s[10] 
                + s[11] * x**s[12] 
                + s[13] * x**s[14] 
                + s[15] * x**s[16]
            )
        ) - (
            (s[1] * s[2] * x**(s[2] - 1))/(x**2 - s[3]**s[4]) 
            - (2 * s[1] * x**(s[2] + 1))/(x**2 - s[3]**s[4])**2 
            + (s[5] * s[6] * x**(s[6] - 1))/(x**2 - s[7]**s[8]) 
            - (2 * s[5] * x**(s[6] + 1))/(x**2 - s[7]**s[8])**2 
            + s[9] * s[10] * x**(s[10] - 1) 
            + s[11] * s[12] * x**(s[12] - 1) 
            + s[13] * s[14] * x**(s[14] - 1) 
            + s[15] * s[16] * x**(s[16] - 1)
        )**2 / (
            4 * (
                s[0] 
                + (s[1] * x**s[2])/(x**2 - s[3]**s[4]) 
                + (s[5] * x**s[6])/(x**2 - s[7]**s[8]) 
                + s[9] * x**s[10] 
                + s[11] * x**s[12] 
                + s[13] * x**s[14] 
                + s[15] * x**s[16]
            )**(3/2)
        )
    return y




def gvd_for_wavelengths(db_record, wl_um, interpolate_kind = "default", verbose = 0):
    """
    db_record is the dictionary
    """
    
    if "type" not in db_record:
        raise KeyError("gvd_for_wavelengths: The database record does not have a key 'type'.")
        
    wl_um = CF.make_numpy_ndarray(wl_um)
    
    error_string = "GVD is not implemented for "
    
    # check the range
    if "formula" in db_record["type"]:   
        if wl_um[0] < db_record["range"][0]: 
            raise ValueError("Error, wavelength %1.2f micron is too low! It should be above %1.2f micron." % (wl_um[0], db_record["range"][0]) )
        elif wl_um[-1] > db_record["range"][1]:
            raise ValueError("Error, wavelength %1.2f micron is too high! It should be below %1.2f micron." % (wl_um[-1], db_record["range"][-1]) )
    
    if "tabulated" in db_record["type"]:   
        if wl_um[0] < db_record["data"][:,0][0]: 
            raise ValueError("Error, wavelength %1.2f micron is too low! It should be above %1.2f micron." % (wl_um[0], db_record["data"][:,0][0]) )
        elif wl_um[-1] > db_record["data"][:,0][-1]:
            raise ValueError("Error, wavelength %1.2f micron is too high! It should be below %1.2f micron." % (wl_um[-1], db_record["data"][:,0][-1]) )
        
    if db_record["type"] == "formula 1":
        DEBUG.verbose("  Using formula 1 to calculate group velocity dispersion", verbose_level = 1)
        gvd = gvd_formula_1(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)

    elif db_record["type"] == "formula 2":
        """
        Formula 2 is the same as formula 1, but some given coefficients are already squared. The square root of these coefficients is taken and the GVD equation for formula 1 is used. 
        """
        DEBUG.verbose("  Using formula 2 to calculate group velocity dispersion", verbose_level = 1)
        s = numpy.copy(db_record["coefficients"])
        for i in range(len(s)):
            if i > 0 and i % 2 == 0:
                s[i] = numpy.sqrt(s[i])
        gvd = gvd_formula_1(wl_um, s)
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)

    elif db_record["type"] == "formula 3":
        DEBUG.verbose("  Using formula 3 to calculate group velocity dispersion", verbose_level = 1)
        gvd = gvd_formula_3(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)
        
    elif db_record["type"] == "formula 4":
        DEBUG.verbose("  Using formula 4 to calculate group velocity dispersion", verbose_level = 1)
        gvd = gvd_formula_4(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)
        
    elif db_record["type"] == "formula 5":
        raise NotImplementedError(error_string + "formula 5")   
        
    elif db_record["type"] == "formula 6":
        raise NotImplementedError(error_string + "formula 6")
        
    elif db_record["type"] == "formula 7":
        raise NotImplementedError(error_string + "formula 7")
        
    elif db_record["type"] == "formula 8":
        raise NotImplementedError(error_string + "formula 8")

    elif db_record["type"] == "formula 9":
        raise NotImplementedError(error_string + "formula 9")

    elif db_record["type"] == "tabulated n":
        raise NotImplementedError(error_string + "tabulated data (type n)")
        
    elif db_record["type"] == "tabulated nk":
        raise NotImplementedError(error_string + "tabulated data (type nk)")
    
    else:
        raise ValueError("The type of data in the database record is unknown (usually formula 1-9 or tabulated data). Type here is %s." % db_record["type"])

    return gvd
        

        
        
        
def get_ri(paf, wl_um = [], um_range = [0.3, 0.6], n_steps = 100, interpolate_kind = "default", verbose = 0):
    """
    
    INPUT:
    paf (str): path and filename
    wl_um (ndarray): wavelength axis, in micrometer
    um_range (list with 2 elements): if wl_um is not given, it will plot a range.
    n_steps (int): number of steps to plot the range.
    ax (plt axis): if False, it will make a new figure
    interpolate_kind (str): for tabulated data, the type of interpolation
    
    """
    temp = paf.split("/")
    DEBUG.verbose("Importing data for %s by %s" % (temp[-2], temp[-1][:-4]), verbose_level = 1)

    db_record = RIRY.import_refractive_index(paf = paf, verbose = verbose)
    DEBUG.verbose("  Imported data", verbose_level = 1)
    
    if len(wl_um) == 0:
        wl_um = numpy.linspace(um_range[0], um_range[1], num = n_steps)
    ri = ri_for_wavelengths(db_record, wl_um, verbose = verbose)

    if type(ri) == int:
        return 0, 0

    return wl_um, ri


def get_gvd(paf, wl_um = [], um_range = [0.3, 0.6], n_steps = 100, interpolate_kind = "default", verbose = 0):
    """
    
    INPUT:
    paf (str): path and filename
    wl_um (ndarray): wavelength axis, in micrometer
    um_range (list with 2 elements): if wl_um is not given, it will plot a range.
    n_steps (int): number of steps to plot the range.
    ax (plt axis): if False, it will make a new figure
    interpolate_kind (str): for tabulated data, the type of interpolation
    
    """
    temp = paf.split("/")
    DEBUG.verbose("Importing data for %s by %s" % (temp[-2], temp[-1][:-4]), verbose_level = 1)

    db_record = RIRY.import_refractive_index(paf = paf, verbose = verbose)
    DEBUG.verbose("  Imported data", verbose_level = 1)

    # error checking
#     del db_record["type"]
#     db_record["type"] = "fiets"

    if len(wl_um) == 0:
        wl_um = numpy.linspace(um_range[0], um_range[1], num = n_steps)
    gvd = gvd_for_wavelengths(db_record, wl_um, verbose = verbose)

    if type(gvd) == int:
        return 0, 0

    return wl_um, gvd


def plot_ri(paf, wl_um = [], um_range = [0.3, 0.6], n_steps = 100, ax = False, interpolate_kind = "default", verbose = 0):
    """
    
    INPUT:
    paf (str): path and filename
    wl_um (ndarray): wavelength axis, in micrometer
    um_range (list with 2 elements): if wl_um is not given, it will plot a range.
    n_steps (int): number of steps to plot the range.
    ax (plt axis): if False, it will make a new figure
    interpolate_kind (str): for tabulated data, the type of interpolation
    
    """
    
    wl_um, ri = get_ri(paf = paf, wl_um = wl_um, um_range = um_range, n_steps = n_steps, interpolate_kind = interpolate_kind, verbose = verbose)
        
    if ax == False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    ax.plot(wl_um, ri)

    ax.set_xlabel("Wavelength (micron)")
    ax.set_ylabel("Index of refraction n")
   
    return ri



def plot_gvd(paf, wl_um = [], um_range = [0.3, 0.6], n_steps = 100, ax = False, interpolate_kind = "default", verbose = 0):
    """
    
    INPUT:
    paf (str): path and filename
    wl_um (ndarray): wavelength axis, in micrometer
    um_range (list with 2 elements): if wl_um is not given, it will plot a range.
    n_steps (int): number of steps to plot the range.
    ax (plt axis): if False, it will make a new figure
    interpolate_kind (str): for tabulated data, the type of interpolation
    
    """
    
    wl_um, gvd = get_gvd(paf, wl_um = wl_um, um_range = um_range, n_steps = n_steps, interpolate_kind = interpolate_kind, verbose = verbose)

    if ax == False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    ax.plot(wl_um, gvd)
    
    ax.set_xlabel("Wavelength (micron)")
    ax.set_ylabel("GVD (fs^2/mm)")
    
    return gvd
      



if __name__ == "__main__": 
    plt.close("all")
    
    path = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndexDB/"
    
#     
#     # formula 1
#     # https://refractiveindex.info/?shelf=main&book=ZnSe&page=Connolly
#     paf = path + "database/main/ZnSe/Connolly.yml"
#     ri = plot_ri(paf, um_range = [4, 18])
#     gvd = plot_gvd(paf, um_range = [4, 8])
#     print(ri[0])
#     print(gvd[0])
    
#     
#     # formula 2
#     # https://refractiveindex.info/?shelf=main&book=CaF2&page=Daimon-20
#     paf = path + "database/main/CaF2/Daimon-20.yml"
#     wl_um = numpy.array([0.3, 0.4, 0.5, 0.6])
# #     plot_ri(paf, wl_um = wl_um)
#     ri = plot_ri(paf)
#     gvd = plot_gvd(paf, um_range = [0.3, 0.6])
#     print(gvd[0])
#     
#     
# #     paf = path + "database/main/CaF2/Daimon-20.yml"
# #     wl_um = numpy.array([0.3, 20])
# #     plot_ri(paf, wl_um = wl_um)
# # #     plot_ri(paf)
#     
#     
#     paf = path + "database/main/Al2O3/Boidin.yml"
# #     wl_um = numpy.array([0.3, 0.4, 0.5, 0.6])
# #     plot_ri(paf, wl_um = wl_um)
#     plot_ri(paf, um_range = [0.4, 17])
# 
# 
#     # tabulated nk
#     paf = path + "database/main/Ag/Babar.yml"
#     plot_ri(paf, um_range = [0.21, 12])
# 
# 
    # formula 3
#     paf = path + "database/organic/C6H6 - benzene/Moutzouris.yml"
#     ri = plot_ri(paf, um_range = [0.5, 1.5])
#     gvd = plot_gvd(paf, um_range = [0.5, 1.5])
#     print(ri[0])
#     print(gvd[0])

    # formula 4
    # https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Eimerl-o
    paf = path + "database/main/BaB2O4/Eimerl-o.yml"
    ri = plot_ri(paf, um_range = [0.25, 1.0])
    gvd = plot_gvd(paf, um_range = [0.25, 1.0])
    print(ri[0])
    print(gvd[0])

    # formula 4
    # https://refractiveindex.info/?shelf=main&book=BiB3O6&page=Umemura-α
    paf = path + "database/main/BiB3O6/Umemura-alpha.yml"
    ri = plot_ri(paf, um_range = [0.5, 1.0])
    gvd = plot_gvd(paf, um_range = [0.5, 1.0])
    print(ri[0])
    print(gvd[0])

 
#     # formula 5
#     # https://refractiveindex.info/?shelf=organic&book=octane&page=Kerl-293K
#     paf = path + "database/organic/C8H18 - octane/Kerl-293K.yml"
#     plot_ri(paf, um_range = [0.33, 0.64])

# 
#     # formula 6
#     # https://refractiveindex.info/?shelf=main&book=H2&page=Peck
#     paf = path + "database/main/H2/Peck.yml"
#     plot_ri(paf, um_range = [0.17, 1.65])







    
    plt.show()
    