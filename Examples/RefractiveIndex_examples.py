"""
Examples on how to use RefractiveIndex.py.

Wavelengths are in microns!

"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import NPCTools.RefractiveIndex as RI

plt.close("all")

# paths. paf = path and filename
path_to_ri_database = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndexDB/"
path_to_file = "database/main/CaF2/Daimon-20.yml"
paf = path_to_ri_database + path_to_file

fig = plt.figure()
n_ax = 2
ax = [0] * n_ax
ax[0] = fig.add_subplot(111)
ax[1] = ax[0].twinx()


# to get specific wavelengths, use wl_um
wl_um = [0.3, 0.6, 0.9]
wl_um, ri = RI.get_ri(paf, wl_um)
wl_um, gvd = RI.get_gvd(paf, wl_um)
for i in range(len(wl_um)):
    print(wl_um[i], ri[i], gvd[i])


# for a plot, you can specify a range
um_range = [0.2, 2.3]
wl_um, ri = RI.get_ri(paf, um_range = um_range)
wl_um, gvd = RI.get_gvd(paf, um_range = um_range)
ax[0].plot(wl_um, ri, color = "blue")
ax[1].plot(wl_um, gvd, color = "red")

ax[0].set_xlabel("Wavelength (micron)")
ax[0].set_ylabel("Index of refraction")
ax[1].set_ylabel("Group velocity dispersion (fs^2/mm)")


plt.show()



if __name__ == "__main__": 
    pass