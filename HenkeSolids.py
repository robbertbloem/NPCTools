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
import NPCTools.Resources.Mathematics as MATH
import NPCTools.Data.HenkeSolidTransmission.HenkeSolidData

imp.reload(NPCTools.Data.HenkeSolidTransmission.HenkeSolidData)

def import_transmission_data(compound_name, verbose = 0):
    """
    Takes a name for a material and outputs the wavelength axis (in eV) and the absorption for 1 micron.
    
    INPUT:
    material: name of material. The data for the gas has to be in '/Data/HenkeSolidTransmission/' and in the python file there. 
    
    OUTPUT:
    - e_ev (ndarray): ev axis
    - tr_norm (ndarray): transmission, normalized to 1 micron
    
    CHANGELOG:
    
    """
    # find path
    path = NPCTools.Data.HenkeSolidTransmission.HenkeSolidData.data_path()
    
    # convert to title case
    compound_name = compound_name.title()
    files_list = NPCTools.Data.HenkeSolidTransmission.HenkeSolidData.files_list()
    
    # find the name in the list
    mat = [item for item in files_list if compound_name in item["compound"]]

    # if not name was found, throw an error
    if len(mat) == 0:
        DEBUG.printError("Compound %s not found" % compound_name, inspect.stack())
        return [0], [0]
    else:
        mat = mat[0]
    
    if verbose > 0:
        s = "HenkeSolids.py:import_transmission_data: Found this item for element %s: %s" % (compound_name, str(mat))
        DEBUG.verbose(s, verbose, 1)
          
    print("Importing data for compound %s" % (compound_name)) 

    # import data
    data = numpy.loadtxt(path + "/" + mat["filename"], skiprows = 2)
    e_ev = data[:,0]
    tr = data[:,1]
    
    # normalize
    ab = numpy.log10(tr)
    ab /= mat["um"]
    tr_norm = 10**ab
    
    return e_ev, tr_norm




def transmission_for_compound(e_ev, tr_norm, ums, evs = [], interpolation_kind = "default"):

    """
    
    INPUTS:
    ev (ndarray): ev axis of the imported data
    
    
    x_ev (ndarray): energies, x-axis
    y_tr (ndarray): absorption, y-axis
    
    ums (list): list with pathlengths, in um
    evs (list): list with specific wavelengths to calculate the transmission for. If length is zero, all energies will be calculated. 
    
    OUTPUT:
    ev (ndarray): 
    
    
    """
    
    ums = CF.make_numpy_ndarray(ums)
    evs = CF.make_numpy_ndarray(evs)

    if len(evs) != 0:
        y_tr = MATH.interpolate_data(original_x = ev, original_y = transmission, new_x = evs, interpolate_kind = interpolation_kind)
        x_ev = evs
        n_evs = len(evs)
    else:
        y_tr = transmission[:]
        x_ev = ev[:]
        n_evs = len(y_tr)

    n_ums = len(ums)
    
    tr = numpy.zeros((n_evs, n_ums))

    for i in range(n_evs):    
        tr[i,:] = y_tr[i]
    
    ab = numpy.log10(tr)
    for i in range(n_ums):    
        ab[:,i] *= ums[i]
    tr = 10**ab
        
    return ev, tr



if __name__ == "__main__": 
    
    
    
    
    interpolation_kind = "linear" # "quadratic" # "linear" # "cubic"
    
    e_ev, tr_norm = import_transmission_data("Zr", verbose = 0)
    ums = [0.1]
    evs = [80, 92]
    e_ev, tr = transmission_for_compound(e_ev, tr_norm, ums, evs = evs, interpolation_kind = "default")
    print(e_ev, tr)
#     i = 8
#     print(ev[i], tr[i])
    
    
#     print(transmission[10])
#     ev, transmission = import_transmission_data("Test", verbose = 0)
#     print(transmission[10])
    
#     mbars = 10
#     cms = 20
#     evs = 30
#     
#     res = absorption_for_gas(ev, absorption, mbars, cms, evs, interpolation_kind = interpolation_kind)
#     
#     print(res)
    
    