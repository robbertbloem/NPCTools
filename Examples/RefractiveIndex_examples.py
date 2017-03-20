"""
Examples on how to use RefractiveIndex.py.

Wavelengths are in microns!
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import imp

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import NPCTools.RefractiveIndex as RI

imp.reload(RI)

plt.close("all")

path_to_ri_database = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndexDB/"

def example_ri_gvd_table():
    """
    Calculate the index of refraction and GVD as a function of wavelengths.
    By using wl_um, the RI and GVD for exactly these wavelengths are returned.
    """
    path_to_file = "database/main/CaF2/Daimon-20.yml"
    paf = path_to_ri_database + path_to_file # paf = path and filename

    # to get specific wavelengths, use wl_um
    wl_um = [0.3, 0.6, 0.9]
    wl_um, ri = RI.get_ri(paf, wl_um)
    wl_um, gvd = RI.get_gvd(paf, wl_um)
    print("WL (um)  RI       GVD (fs^2/mm)")
    for i in range(len(wl_um)):
        print("%2.1f      %5.4f   %5.4f" % (wl_um[i], ri[i], gvd[i]))

 
 
def example_ri_gvd_plot():
    """
    Calculate the index of refraction and GVD as a function of wavelengths.
    um_range gives the start and end of the range to be calculated. 
    """
    path_to_file = "database/main/CaF2/Daimon-20.yml"
    paf = path_to_ri_database + path_to_file # paf = path and filename

    # initialize plot
    fig = plt.figure()
    n_ax = 2
    ax = [0] * n_ax
    ax[0] = fig.add_subplot(211)
    ax[1] = fig.add_subplot(212) 

    # for a plot, you can specify a range
    um_range = [0.2, 2.3]
    
    # RI  
    ax_i = 0
    wl_um, ri = RI.get_ri(paf, um_range = um_range)
    ax[ax_i].plot(wl_um, ri, color = "blue")
    ax[ax_i].set_ylabel("Index of refraction")
    
    # GVD
    ax_i = 1
    wl_um, gvd = RI.get_gvd(paf, um_range = um_range)
    ax[ax_i].plot(wl_um, gvd, color = "red")
    ax[ax_i].set_ylabel("GVD (fs^2/mm)")
    
    ax[0].set_title("Index of refraction and GVD for CaF2")
    ax[1].set_xlabel("Wavelength (micron)")   
    

def effect_thickness_on_gvd():
    """
    Compare the effects of two different thicknesses on the pulse length.
    """
    path_to_file = "database/glass/schott/N-BK7.yml"
    paf = path_to_ri_database + path_to_file

    # initialize plot
    fig = plt.figure()
    n_ax = 1
    ax = [0] * n_ax
    ax[0] = fig.add_subplot(111)

    # pulse length as a function of wavelength, for 2 thicknesses
    um_range = [0.3, 2.5]
    wl_um, gvd = RI.get_gvd(paf, um_range = um_range)
    # t_fs and d_mm can be a number, a list or an ndarray
    t_fs = 35 
    d_mm = [5, 10]
    # t_out is a 3D matrix with shape wavelengths x (initial) t_fs x d_mm
    t_out = RI.effect_gvd(wl_um, gvd, t_fs, d_mm)
    ax_i = 0
    ax[ax_i].plot(wl_um, t_out[:,0,0], label = "5 mm")
    ax[ax_i].plot(wl_um, t_out[:,0,1], label = "10 mm")
    ax[ax_i].legend()
    ax[ax_i].set_title("Initial pulse length is 35 fs, for BK7 glass")
    ax[ax_i].set_xlabel("Wavelength (micron)")
    ax[ax_i].set_ylabel("Pulse length (fs)")


def effect_gvd_for_pulselengths():
    """
    Plot the pulse length before and after a material. This is done for two materials.
    """
    # initialize plot
    fig = plt.figure()
    n_ax = 1
    ax = [0] * n_ax
    ax[0] = fig.add_subplot(111)
    ax_i = 0
    
    # use only 1 wavelength
    wl_um = [0.8]
    # make an array of initial times. ndarray: [10, 11, 12, ..., 149]
    t_fs = numpy.arange(5,100)
    # only one thickness
    d_mm = [10]
    
    # BK7
    path_to_file = "database/glass/schott/N-BK7.yml"
    paf = path_to_ri_database + path_to_file
    wl_um, gvd = RI.get_gvd(paf, wl_um = wl_um) 
    t_out = RI.effect_gvd(wl_um, gvd, t_fs, d_mm)
    ax[ax_i].plot(t_fs, t_out[0,:,0], label = "BK7")
    
    # CaF2
    path_to_file = "database/main/CaF2/Daimon-20.yml"
    paf = path_to_ri_database + path_to_file
    wl_um, gvd = RI.get_gvd(paf, wl_um = wl_um) 
    t_out = RI.effect_gvd(wl_um, gvd, t_fs, d_mm)
    ax[ax_i].plot(t_fs, t_out[0,:,0], label = "CaF2")    
    
    # plot
    ax[ax_i].legend()
    ax[ax_i].grid()
    ax[ax_i].set_xlim(0,100)
    ax[ax_i].set_ylim(0,125)
    # diagonal line
    ax[ax_i].plot([0,200], [0,200], color = "black", linewidth = 1)
    ax[ax_i].plot([35,35], [0,200], color = "grey", linewidth = 1)
    ax[ax_i].set_title("Effect on pulse length for 10 mm material at 800 nm")
    ax[ax_i].set_xlabel("Initial pulse length (fs)")
    ax[ax_i].set_ylabel("Resulting pulse length (fs)")
    
    





if __name__ == "__main__": 
    example_ri_gvd_table()
    example_ri_gvd_plot()
    effect_thickness_on_gvd()
    effect_gvd_for_pulselengths()
    
    plt.show()