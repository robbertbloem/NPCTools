from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect
import re

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import NPCTools.Debug as DEBUG

def make_numpy_ndarray(val):
    """
    Make a numpy.ndarray out of val. 
    
    Types of val that are accepted:        
    int, float, string: make it a list and then numpy.ndarray.
    list: make it a numpy.ndarray
    tuple: convert to a list, then numpy.ndarray
    numpy.ndarray: return directly
    
    Not accepted:
    dict
    
    CHANGELOG:
    20130317/RB: started
    20170310/RB: raise a type error for dict, instead of my own error
    
    """
            
    if type(val) == numpy.ndarray:
        return val
    elif type(val) == list:
        return numpy.array(val)
    elif type(val) == dict:
        raise TypeError("Value can't be a dict")
#         DEBUG.printError("Value can't be a dict or tuple", inspect.stack())
#         return False
    elif type(val) == tuple:
        return numpy.array(list(val))
    else:
        return numpy.array([val])  


def string_with_numbers_to_list(string):
    """
    Receive a string with numbers, for example from a file, and make a ndarray out of it. The output type is always float. Commas and spaces indicate the separation between numbers and can be mixed. Newlines are removed.
    
    string = "0, 0.1,\n 1e+3 1e-2"
    output: [0.0, 0.1, 1000, 0.01]

    CHANGELOG:
    20170309/RB: started

    """
    
    string = string.replace("\n", " ")
    # remove everything, except \d (numbers), \s (spaces), . (decimal), e (exponent), +, - (signs)
    non_decimal = re.compile(r'[^\d\s.e+-]+')
    res = non_decimal.sub('', string)
    data = res.split(" ")
    # remove excess spaces
    data = list(filter(None, data))

    data = numpy.array(data)
    data = data.astype(float)

    return data
    
    







if __name__ == "__main__": 
    pass