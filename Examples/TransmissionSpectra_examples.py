from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import imp

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import NPCTools.TransmissionSpectra as TS

imp.reload(TS)

plt.close("all")


def simple_example():
    
    compound_name = "sapphire"
    
    wl_nm, abs_mm, Rs, Rp = TS.import_transmission_data(compound_name, verbose = 0)
    thickness_mms = [10, 20]
    wl_nms, thickness_mms, transmission = TS.transmission_for_compound(wl_nm, abs_mm, R = Rs, thickness_mms = thickness_mms, wl_nms = [], interpolation_kind = "default")
    
    # initialize plot
    fig = plt.figure()
    n_ax = 1
    ax = [0] * n_ax
    ax[0] = fig.add_subplot(111)
    
    ax_i = 0
    ax[ax_i].plot(wl_nms, transmission[0,:])
    ax[ax_i].plot(wl_nms, transmission[1,:])



def absorption_transmission_example():
    
    compound_name = "uvfs"
    
    wl_nm, abs_mm, Rs, Rp = TS.import_transmission_data(compound_name, wl_range = (0.21, 3.71), verbose = 0)
    thickness_mms = [10]
    wl_nms, thickness_mms, transmission = TS.transmission_for_compound(wl_nm, abs_mm, R = Rs, thickness_mms = thickness_mms, wl_nms = [], interpolation_kind = "default")
    
    # initialize plot
    fig = plt.figure()
    n_ax = 1
    ax = [0] * n_ax
    ax[0] = fig.add_subplot(111)

    
    ax_i = 0
    ax[ax_i].plot(wl_nms, transmission[0,:], label = "Transmission")
    ax[ax_i].plot(wl_nms, transmission[0,:] / (1 - Rs)**2, label = "Absorption component")
    ax[ax_i].plot(wl_nms, (1 - Rs)**2, label = "Reflection component")

    ax[ax_i].legend()



if __name__ == "__main__": 
#     simple_example()
    absorption_transmission_example()
    
    
    plt.show()