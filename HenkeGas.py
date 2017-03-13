from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import imp
import os
import inspect

import numpy
import matplotlib 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import NPCTools.Debug as DEBUG
import NPCTools.Resources.CommonFunctions as CF
import NPCTools.Data.HenkeGasTransmission.HenkeGasData

imp.reload(NPCTools.Data.HenkeGasTransmission.HenkeGasData)

def import_absorption_data(gas_name, verbose = 0):
    """
    Takes a name for a gas and outputs the wavelength axis (in eV) and the absorption for 1 mbar and 1 cm path length.
    
    INPUT:
    gas_name: name of gas. The data for the gas has to be in '/Data/HenkeGasTransmission/' and in the python file there. 
    
    OUTPUT:
    - ev (ndarray): ev axis
    - absorption (ndarray): absorption, normalized to 1 mbar and 1 cm
    
    CHANGELOG:
    
    """
    # find path
    path = NPCTools.Data.HenkeGasTransmission.HenkeGasData.data_path()
    
    # convert to title case
    gas_name = gas_name.title()
    files_list = NPCTools.Data.HenkeGasTransmission.HenkeGasData.files_list()
    gas = [item for item in files_list if item["gas"] == gas_name]

    if len(gas) == 0:
        DEBUG.printError("Gas %s not found" % gas_name, inspect.stack())
        return [0], [0]
    else:
        gas = gas[0]
        
    if verbose > 0:
        s = "HenkeGas.py:import_absorption_data: Found this item for gas %s: %s" % (gas_name, str(gas))
        DEBUG.verbose(s, verbose, 1)
          
    print("Importing data for %s gas" % (gas_name)) 

    data = numpy.loadtxt(path + "/" + gas["filename"], skiprows = 2)
    ev = data[:,0]
    transmission = data[:,1]
    absorption = numpy.log10(transmission)

    # normalize to 1 mbar and 1 cm
    absorption /= gas["mbar"]
    absorption /= gas["cm"]
    
    return ev, absorption


def interpolate_absorption(ev, ab, x, interpolation_kind = ""):
    """
    filter_data is the data
    x is the value we want
    """    
    
    if interpolation_kind == "":
        interpolation_kind = "linear"

    f = interp1d(ev, ab, kind = interpolation_kind)
    absorption = f(x)    

    return absorption


def absorption_for_gas(ev, absorption, mbars, cms, evs = [], interpolation_kind = ""):

    """
    
    INPUTS:
    x_ev (ndarray): energies, x-axis
    y_ab (ndarray): absorption, y-axis
    
    mbars (list): list with pressures, in mbar
    cms (list): list with pathlengths, in cm
    evs (list): list with specific wavelengths to calculate the transmission for. If length is zero, all energies will be calculated. 
    
    OUTPUT:
    ev (ndarray): 
    
    
    """
    
    mbars = CF.make_numpy_ndarray(mbars)
    cms = CF.make_numpy_ndarray(cms)
    evs = CF.make_numpy_ndarray(evs)

    if len(evs) != 0:
        y_ab = interpolate_absorption(ev = ev, ab = absorption, x = evs, interpolation_kind = interpolation_kind)
        x_ev = evs
        n_evs = len(evs)
    else:
        y_ab = absorption[:]
        x_ev = ev[:]
        n_evs = len(y_ab)
    
    n_mbars = len(mbars)
    n_cms = len(cms)
    
    ab = numpy.zeros((n_evs, n_mbars, n_cms))

    for i in range(n_evs):    
        ab[i,:,:] = y_ab[i]
    
    for i in range(n_mbars):    
        ab[:,i,:] *= mbars[i]
    
    for i in range(n_cms):    
        ab[:,:,i] *= cms[i]
        
    tr = 10**ab
        
    return ev, ab, tr



if __name__ == "__main__": 
    
    interpolation_kind = "linear" # "quadratic" # "linear" # "cubic"
    
    ev, absorption = import_absorption_data("helium", verbose = 1)
    
    mbars = 10
    cms = 20
    evs = 30
    
    res = absorption_for_gas(ev, absorption, mbars, cms, evs, interpolation_kind = interpolation_kind)
    
    print(res)
    
    