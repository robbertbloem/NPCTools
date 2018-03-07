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
import NPCTools.Data.TransmissionSpectra.TransmissionSpectraData

import NPCTools.RefractiveIndex as RI

imp.reload(RI)

imp.reload(NPCTools.Data.TransmissionSpectra.TransmissionSpectraData)

def import_transmission_data(compound_name, wl_range = (), verbose = 0):
    """
    Takes a name for a material and outputs the wavelength axis (in eV) and the absorption for 1 micron.
    
    INPUT:
    material: name of material. The data for the gas has to be in '/Data/TransmissionSpectra/' and in the python file there. 
    
    OUTPUT:
    - wl_nm (ndarray): wavelength axis
    - tr_norm (ndarray): transmission, normalized to 1 mm
    
    CHANGELOG:
    
    """
    # find path
    path = NPCTools.Data.TransmissionSpectra.TransmissionSpectraData.data_path()
    
    # convert to title case
    compound_name = compound_name.lower()
    files_list = NPCTools.Data.TransmissionSpectra.TransmissionSpectraData.files_list()
    
    # find the name in the list
    mat = [item for item in files_list if compound_name in item["compound"]]

    # if not name was found, throw an error
    if len(mat) == 0:
        DEBUG.printError("Compound %s not found" % compound_name, inspect.stack())
        return [0], [0]
    else:
        mat = mat[0]
    
    if verbose > 0:
        s = "TransmissionSpectra.py:import_transmission_data: Found this item for element %s: %s" % (compound_name, str(mat))
        DEBUG.verbose(s, verbose, 1)
          
    print("Importing data for compound %s" % (compound_name)) 

    # import data
    data = numpy.loadtxt(path + "/" + mat["filename"], comments = "#", delimiter = ",")
    wl = data[:,0]
    tr = data[:,1]
    
    if wl[0] > wl[-1]:
        wl = wl[::-1]
        tr = tr[::-1]
    
    if mat["wl_unit"] == "nm":
        wl_nm = wl
        wl_um = wl / 1000
    elif mat["wl_unit"] == "um":
        wl_nm = 1000 * wl
        wl_um = wl

    temp = numpy.where(wl_um <= wl_range[0])[0]
    if len(temp) > 0:
        wl_nm = wl_nm[temp[-1]:]
        wl_um = wl_um[temp[-1]:]
        tr = tr[temp[-1]:]

    temp = numpy.where(wl_um >= wl_range[1])[0]
    if len(temp) > 0:
        wl_nm = wl_nm[:temp[0]]
        wl_um = wl_um[:temp[0]]
        tr = tr[:temp[0]]
    
    a_deg = [0]
    
    path_to_ri_database = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndexDB/"
    paf = path_to_ri_database + mat["refractive_index_path"]
    
    # air -> material
    wl_um, a_deg, Rs, Rp = RI.get_reflectance(paf2 = paf, wl_um = wl_um, a_deg = a_deg)

    Rs = Rs[0,:]
    Rp = Rp[0,:]

    if mat["tr_unit"] == "pct":
        tr /= 100

    
    R = (1-Rs)**2
    absorption_component = tr / R #(1 - R)
    abs = -numpy.log10(absorption_component)
    abs_mm = abs / mat["mm"]
    
    numpy.putmask(abs_mm, abs_mm < 0, 0)

    
    return wl_nm, abs_mm, Rs, Rp

def import_transmission_data_raw(compound_name, verbose = 0):
    """
    Takes a name for a material and outputs the wavelength axis (in eV) and the absorption for 1 micron.
    
    INPUT:
    material: name of material. The data for the gas has to be in '/Data/TransmissionSpectra/' and in the python file there. 
    
    OUTPUT:
    - wl_nm (ndarray): wavelength axis
    - tr_norm (ndarray): transmission, normalized to 1 mm
    
    CHANGELOG:
    
    """
    # find path
    path = NPCTools.Data.TransmissionSpectra.TransmissionSpectraData.data_path()
    
    # convert to title case
    compound_name = compound_name.lower()
    files_list = NPCTools.Data.TransmissionSpectra.TransmissionSpectraData.files_list()
    
    # find the name in the list
    mat = [item for item in files_list if compound_name in item["compound"]]

    # if not name was found, throw an error
    if len(mat) == 0:
        DEBUG.printError("Compound %s not found" % compound_name, inspect.stack())
        return [0], [0]
    else:
        mat = mat[0]
    
    if verbose > 0:
        s = "TransmissionSpectra.py:import_transmission_data: Found this item for element %s: %s" % (compound_name, str(mat))
        DEBUG.verbose(s, verbose, 1)
          
    print("Importing data for compound %s" % (compound_name)) 

    # import data
    data = numpy.loadtxt(path + "/" + mat["filename"], comments = "#", delimiter = ",")
    wl = data[:,0]
    tr = data[:,1]
    
    if mat["wl_unit"] == "nm":
        wl_nm = wl
        wl_um = wl / 1000
    elif mat["wl_unit"] == "um":
        wl_nm = 1000 * wl
        wl_um = wl
    
    if mat["tr_unit"] == "pct":
        tr /= 100

    return wl_nm, tr


def transmission_for_compound(wl_nm, abs_mm, R, thickness_mms, wl_nms = [], interpolation_kind = "default"):

    """
    
    INPUTS:
    wl_nm (ndarray): wavelength axis of the imported data
    abs_mm (ndarray): pure absorption component
    R (ndarray): reflection component
    
    thickness_mms (int, list, ndarray): list with thicknesses in mm
    wl_nms (int, list, ndarray): list with wavelengths of interest, can be empty, then wl_nm is used. 

    mms (list): list with pathlengths, in mm
    l_nms (list): list with specific wavelengths to calculate the transmission for. If length is zero, all energies will be calculated. 
    
    OUTPUT:
    ev (ndarray): 
    
    
    """
    
    thickness_mms = CF.make_numpy_ndarray(thickness_mms)
    n_thick = len(thickness_mms)
    wl_nms = CF.make_numpy_ndarray(wl_nms)

    if len(wl_nms) != 0:
        y_abs_mm = MATH.interpolate_data(original_x = wl_nm, original_y = abs_mm, new_x = wl_nms, interpolate_kind = interpolation_kind)
        R = MATH.interpolate_data(original_x = wl_nm, original_y = R, new_x = wl_nms, interpolate_kind = interpolation_kind)
        n_nms = len(wl_nms)
    else:
        y_abs_mm = abs_mm[:]
        wl_nms = wl_nm[:]
        n_nms = len(abs_mm)

    ABS, TH = numpy.meshgrid(y_abs_mm, thickness_mms)
    R, dump = numpy.meshgrid(R, thickness_mms)
    
    _R = (1-R)**2    
    transmission = 10**(-ABS * TH) * (_R)

        
    return wl_nms, thickness_mms, transmission








if __name__ == "__main__": 
    
    plt.close("all")
    
    
    interpolation_kind = "linear" # "quadratic" # "linear" # "cubic"
    
    wl_nm, abs_mm, Rs, Rp = import_transmission_data("sapphire", verbose = 0)
    mms = [1,10]
    wl_nms = [] 
    wl_nms, thickness_mms, transmission = transmission_for_compound(wl_nm, abs_mm = abs_mm, R = Rs, thickness_mms = mms, wl_nms = wl_nms, interpolation_kind = "default")
    print(wl_nms, transmission)

    plt.plot(wl_nms, transmission[0,:])
    plt.plot(wl_nms, transmission[1,:])
    plt.show()
    