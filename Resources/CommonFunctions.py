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
    non_decimal = re.compile(r'[^\d\s.eE+-]+')
    res = non_decimal.sub('', string)
    data = res.split(" ")
    # remove excess spaces
    data = list(filter(None, data))

    data = numpy.array(data)
    data = data.astype(float)

    return data
    
    

def find_index_for_value(list, value, round = "norm"):
    """
    Find the index for a value in a list.
    
    round: norm, up or down
    
    Take a list [10, 11, 12, 13] and a value 11.6. The index can be 1 or 2. 
    - round norm: 12 is closer than 11, so returned index is of 12 (i.e. 2)
    - round up: return the index of 12 (i.e. 2)
    - round down: return the index of 11 (i.e. 1)
     
    """


    list = numpy.array(list)
    idx = numpy.argwhere(list > value)[0][0]
    
    if idx == 0:
        return idx

    v0 = list[idx - 1]
    v1 = list[idx]
    
    if v0 == value:
        return idx - 1
    
    if round == "norm":
        if abs(v0 - value) < abs(v1 - value):
            return idx - 1
        else:
            return idx

    elif round == "up":
        return idx    

    elif round == "down":
        return idx - 1


def find_indices_for_range(list, range, frame = "data", round = "maximize"):

    if frame == "data":
        if round == "maximize": 
            round_a = "down"
            round_b = "up"
        else:
            round_a = "norm"
            round_b = "norm"
    
        a = find_index_for_value(list, range[0], round = round_a)
        b = find_index_for_value(list, range[1], round = round_b)
    else:
        a = range[0]
        b = range[1]
            
    return a, b




if __name__ == "__main__": 
    pass