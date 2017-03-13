from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import os

import numpy
import matplotlib 
import matplotlib.pyplot as plt

def data_path():
    return os.path.split(os.path.realpath(__file__))[0]
#     return __file__
    
    

def files_list():
    l = [
        {"filename": "Ar_1mbar_1cm_295C.dat", "mbar": 1.0, "cm": 1.0, "gas": "Argon"}, 
        {"filename": "Ne_10mbar_1cm_295C.dat", "mbar": 10.0, "cm": 1.0, "gas": "Neon"}, 
        {"filename": "He_10mbar_1cm_295C.dat", "mbar": 10.0, "cm": 1.0, "gas": "Helium"}, 
    ]
    
    return l

if __name__ == "__main__": 
    pass

#     print(data_path())

