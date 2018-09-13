"""

"""

import numpy
import matplotlib 
import matplotlib.pyplot as plt


def gaussian(t, mu, sigma, amplitude, offset):
    return amplitude * numpy.exp( -(t - mu)**2 / (2 * sigma**2) ) + offset




class Fitting():
    
    def __init__(self):
        self.function = None
        self.xdata = []
        self.ydata = []

    def set_xdata(self, xdata):
        self.xdata = xdata

    def set_ydata(self, ydata):
        self.ydata = ydata
        
    def set_function(self, function):
        self.function = function
    
    def calculate_from_params(self, params):
        return self.function(self.xdata, *params)
    





if __name__ == "__main__": 
    pass