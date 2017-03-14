from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import argparse
import unittest
import inspect

from imp import reload

import numpy
import matplotlib 
import matplotlib.pyplot as plt


import NPCTools.Resources.Mathematics as MATH
import NPCTools.Debug as DEBUG


reload(MATH)
reload(DEBUG)





class Test_interpolate_data(unittest.TestCase):
    """
    Some basic checks of the interpolation routine. It checks if it works, but doesn't check the data. 
    """


    def setUp(self):
        self.verbose = 0

    def test_default(self):
        original_x = [0,1]
        original_y = [2,3]
        new_x = numpy.array([0, 0.5, 1.0])
        ri = MATH.interpolate_data(original_x, original_y, new_x, interpolate_kind = "default", verbose = self.verbose)
        ri_check = numpy.array([2.0, 2.5, 3.0])
        self.assertTrue(numpy.all(ri == ri_check))
        
    def test_linear_1(self):     
        original_x = [0,1]
        original_y = [2,3]
        new_x = numpy.array([0, 0.5, 1.0])
        ri = MATH.interpolate_data(original_x, original_y, new_x, interpolate_kind = "default", verbose = self.verbose)
        ri_check = numpy.array([2.0, 2.5, 3.0])
        self.assertTrue(numpy.all(ri == ri_check))
        
    def test_linear_2(self):

        n = 101
        data = numpy.zeros((n, 2))   
        original_x = numpy.arange(n) / 15
        original_y = numpy.sin(original_x)  
        cut_x = original_x[::5]
        cut_y = original_y[::5]        
        ri = MATH.interpolate_data(cut_x, cut_y, original_x, interpolate_kind = "linear", verbose = self.verbose)


    def test_cubic_1(self):

        n = 101
        data = numpy.zeros((n, 2))   
        original_x = numpy.arange(n) / 15
        original_y = numpy.sin(original_x)  
        cut_x = original_x[::5]
        cut_y = original_y[::5]        
        ri = MATH.interpolate_data(cut_x, cut_y, original_x, interpolate_kind = "cubic", verbose = self.verbose)
        
    def test_quadratic_1(self):

        n = 101
        data = numpy.zeros((n, 2))   
        original_x = numpy.arange(n) / 15
        original_y = numpy.sin(original_x)  
        cut_x = original_x[::5]
        cut_y = original_y[::5]        
        ri = MATH.interpolate_data(cut_x, cut_y, original_x, interpolate_kind = "quadratic", verbose = self.verbose)
        

if __name__ == '__main__': 
    

 
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_interpolate_data)
    unittest.TextTestRunner(verbosity=1).run(suite) 
 
 
 



