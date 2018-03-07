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
    wl_um = [0.2]
    # make an array of initial times. ndarray: [10, 11, 12, ..., 149]
    t_fs = numpy.arange(5,100)
    # only one thickness
    d_mm = [0.5]
#     
#     # BK7
#     path_to_file = "database/glass/schott/N-BK7.yml"
#     paf = path_to_ri_database + path_to_file
#     wl_um, gvd = RI.get_gvd(paf, wl_um = wl_um) 
#     t_out = RI.effect_gvd(wl_um, gvd, t_fs, d_mm)
#     ax[ax_i].plot(t_fs, t_out[0,:,0], label = "BK7")
#     
#     # CaF2
#     path_to_file = "database/main/CaF2/Daimon-20.yml"
#     paf = path_to_ri_database + path_to_file
#     wl_um, gvd = RI.get_gvd(paf, wl_um = wl_um) 
#     t_out = RI.effect_gvd(wl_um, gvd, t_fs, d_mm)
#     ax[ax_i].plot(t_fs, t_out[0,:,0], label = "CaF2")    
    
    # Quartz
    d_mm = [0.2 * numpy.sqrt(2)]
    path_to_file = "database/main/SiO2/Ghosh-e.yml"
    paf = path_to_ri_database + path_to_file
    wl_um, gvd = RI.get_gvd(paf, wl_um = wl_um) 
    t_out = RI.effect_gvd(wl_um, gvd, t_fs, d_mm)
    ax[ax_i].plot(t_fs, t_out[0,:,0], label = "Quartz, 0.28 mm")        

    # Sapphire
    d_mm = [0.5]
    path_to_file = "database/main/Al2O3/Malitson-e.yml"
    paf = path_to_ri_database + path_to_file
    wl_um, gvd = RI.get_gvd(paf, wl_um = wl_um) 
    t_out = RI.effect_gvd(wl_um, gvd, t_fs, d_mm)
    ax[ax_i].plot(t_fs, t_out[0,:,0], label = "Sapphire, 0.5 mm") 

    # Sapphire
    d_mm = [0.2]
    path_to_file = "database/main/Al2O3/Malitson-e.yml"
    paf = path_to_ri_database + path_to_file
    wl_um, gvd = RI.get_gvd(paf, wl_um = wl_um) 
    t_out = RI.effect_gvd(wl_um, gvd, t_fs, d_mm)
    ax[ax_i].plot(t_fs, t_out[0,:,0], label = "Sapphire, 0.2 mm") 

    # plot
    ax[ax_i].legend()
    ax[ax_i].grid()
    ax[ax_i].set_xlim(0,60)
    ax[ax_i].set_ylim(0,80)
    # diagonal line
    ax[ax_i].plot([0,200], [0,200], color = "black", linewidth = 1)
    ax[ax_i].plot([35,35], [0,200], color = "grey", linewidth = 1)
#     ax[ax_i].set_title("Effect on pulse length for {:2.1f} mm material at 800 nm".format(d_mm[0]))
    ax[ax_i].set_title("Effect on pulse length at {:3d} nm".format(int(1000 * wl_um[0])))
    ax[ax_i].set_xlabel("Initial pulse length (fs)")
    ax[ax_i].set_ylabel("Resulting pulse length (fs)")
    
    

def effect_reflection_for_wavelengths():
    """
    Compare the effects of two different thicknesses on the pulse length.
    """
    path_to_file = "database/main/CaF2/Daimon-20.yml"
    paf = path_to_ri_database + path_to_file

    # initialize plot
    fig = plt.figure()
    n_ax = 1
    ax = [0] * n_ax
    ax[0] = fig.add_subplot(111)

    # wavelength range
    um_range = [0.21, 2.3]

    a_deg = [0]
    # air -> caf2
    wl_um, a_deg, Rs1, Rp1 = RI.get_reflectance(paf2 = paf, um_range = um_range, a_deg = a_deg)
    # caf2 -> air
    wl_um, a_deg, Rs2, Rp2 = RI.get_reflectance(paf1 = paf, um_range = um_range, a_deg = a_deg)
    Rs = Rs1 + Rs2

    ax_i = 0
    _a = 0

    ax[ax_i].plot(wl_um, 100 * Rs1[_a,:], label = "air -> CaF2")
    ax[ax_i].plot(wl_um, 100 * Rs2[_a,:], label = "CaF2 -> vacuum")
    ax[ax_i].plot(wl_um, 100 * Rs[_a,:], label = "total")

    
    ax[ax_i].legend()
    ax[ax_i].set_title("Reflection for 0 deg AOI for CaF2")
    ax[ax_i].set_xlabel("Wavelength (micron)")
    ax[ax_i].set_ylabel("Reflection (%)")
    ax[ax_i].set_ylim(0,8)


def effect_reflection_for_angles():

    path_to_file = "database/main/CaF2/Daimon-20.yml"
    paf = path_to_ri_database + path_to_file

    # initialize plot
    fig = plt.figure()
    n_ax = 2
    ax = [0] * n_ax
    ax[0] = fig.add_subplot(121)
    ax[1] = fig.add_subplot(122)

    # wavelength range
    wl_um = [1]

    a_range = (0, 90)
    a_deg = numpy.linspace(0, 90, 1000)
    # air -> caf2
    wl_um, a_deg, Rs1, Rp1 = RI.get_reflectance(paf2 = paf, wl_um = wl_um, a_deg = a_deg) #a_range = a_range)
    # caf2 -> air
    wl_um, a_deg, Rs2, Rp2 = RI.get_reflectance(paf1 = paf, wl_um = wl_um, a_deg = a_deg) #a_range = a_range) #a_deg = a_deg)


    
    _wl = 0

    ax_i = 0
    ax[ax_i].plot(a_deg, 100 * Rs1[:,_wl], label = "s-pol")
    ax[ax_i].plot(a_deg, 100 * Rp1[:,_wl], label = "p-pol")
    ax[ax_i].set_title("air -> CaF2")
    
    
    
    ax_i = 1
    ax[ax_i].plot(a_deg, 100 * Rs2[:,_wl], label = "s-pol")
    ax[ax_i].plot(a_deg, 100 * Rp2[:,_wl], label = "p-pol")
    ax[ax_i].set_title("CaF2 -> air")

    for ax_i in range(2):
        ax[ax_i].legend()
        ax[ax_i].set_ylabel("Reflection (%)")
        ax[ax_i].set_xlabel("Angle (degrees)")
        
        ax[ax_i].set_xlim(0,90)
        ax[ax_i].set_ylim(0,100)
        
    
    fig.suptitle("Reflection for 1 micron light for CaF2")




if __name__ == "__main__": 
#     example_ri_gvd_table()
#     example_ri_gvd_plot()
#     effect_thickness_on_gvd()
    effect_gvd_for_pulselengths()
#     effect_reflection_for_wavelengths()
#     effect_reflection_for_angles()
    
    plt.show()