from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import os

import numpy
import matplotlib 
import matplotlib.pyplot as plt

def data_path():
    return os.path.split(os.path.realpath(__file__))[0]

    
    

def files_list():
    l = [
        {"filename": "Zr_100nm.dat", "um": 0.1, "compound": ["Zr", "Zirconium"]}, 
        {"filename": "Zr_300nm.dat", "um": 0.3, "compound": ["Test"]}, 
    ]
    
    return l

if __name__ == "__main__": 
    pass



