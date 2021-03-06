from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from imp import reload

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.PlottingTools as PT
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

global suite_list
suite_list = [
    "Suite 1: Make coordinates",
    "Suite 2: Find longest list",
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")
parser.add_argument("-s1", "--skip1", action = "store_true", help = suite_list[0])
parser.add_argument("-s2", "--skip2", action = "store_true", help = suite_list[1])

# process
args = parser.parse_args()

# reload
if args.reload:
    reload(DEBUG)
    reload(PT)

def execute(args):

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_make_coordinates)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    else:
        DEBUG.verbose("Skipping: " + suite_list[0], True)

    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_longest_list)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    else:
        DEBUG.verbose("Skipping: " + suite_list[1], True)


class Test_make_coordinates(unittest.TestCase):
    """
    See if making coordinates works

    CHANGELOG:
    20130317/RB: started the suite

    """

    def setUp(self):
        """
        Toggle verbose using the command line. 
        """

        self.flag_verbose = args.verbose

    def test_single_plot_10(self):
        """
        Simple test
        Single plot in middle
        inch_per_unit = 1
        """
        
        inch_per_unit = 1
        x_units = 4
        y_units = 4
        left = 1
        bottom = 1
        width = 2
        height = 2
        
        figsize, coords = PT.make_coordinates(inch_per_unit, x_units, y_units, left, bottom, width, height)
        
        self.assertEqual(figsize, (4,4))
        self.assertEqual(coords, [(0.25, 0.25, 0.5, 0.5)])
        

    def test_single_plot_05(self):
        """
        See if changing inch_per_init works
        Single plot in middle
        inch_per_unit = 0.5
        """
        inch_per_unit = 0.5
        x_units = 4
        y_units = 4
        left = 1
        bottom = 1
        width = 2
        height = 2
        
        figsize, coords = PT.make_coordinates(inch_per_unit, x_units, y_units, left, bottom, width, height)
        
        self.assertEqual(figsize, (2,2))
        self.assertEqual(coords, [(0.25, 0.25, 0.5, 0.5)])
    

    def test_double_h_plot_1(self):
        """
        Two plots, equal sizes
        |1<2>2<2>1| is |1/8<1/4>1/4<1/4>1/8|
        """
        inch_per_unit = 1.0
        x_units = 8
        y_units = 4
        left = [1,5]
        bottom = 1
        width = 2
        height = 2
        
        figsize, coords = PT.make_coordinates(inch_per_unit, x_units, y_units, left, bottom, width, height)
                
        self.assertEqual(figsize, (8,4))
        self.assertEqual(coords, [(0.125, 0.25, 0.25, 0.5), (0.625, 0.25, 0.25, 0.5)])




    def test_double_h_plot_2(self):
        """
        Two plots
        |1<3>1<4>1| is |0.1<0.3>0.1<0.4>0.1|

        """
        inch_per_unit = 1.0
        x_units = 10
        y_units = 4
        left = [1,5]
        bottom = 1
        width = [3,4]
        height = 2
        
        figsize, coords = PT.make_coordinates(inch_per_unit, x_units, y_units, left, bottom, width, height)
                
        self.assertEqual(figsize, (10,4))
        self.assertEqual(coords, [(0.1, 0.25, 0.3, 0.5), (0.5, 0.25, 0.4, 0.5)])


    def test_four_plots_simple(self):
        """
        Four plots - simple
        01234
        3 x x
        2xxxx
        1 x x 
        0xxxx
    
        """
        inch_per_unit = 1.0
        x_units = 5
        y_units = 5
        left = [1,3]
        bottom = [3,3,1,1]
        width = 1
        height = 1
        
        figsize, coords = PT.make_coordinates(inch_per_unit, x_units, y_units, left, bottom, width, height)
    
        c = [(0.2, 0.6, 0.2, 0.2), (0.6, 0.6, 0.2, 0.2), (0.2, 0.2, 0.2, 0.2), (0.6, 0.2, 0.2, 0.2)]
        
        self.assertEqual(figsize, (5,5))
        self.assertTrue(numpy.allclose(numpy.array(coords), numpy.array(c)))


    def test_four_plots_complex(self):
        """
        Four plots - quite complex
        01234567
        6   x  x
        5   x  x
        4xxxx  x
        3  xxxxx
        2  x   x
        1  x   x
        0xxxxxxx
    
        """
        inch_per_unit = 1.0
        x_units = 8
        y_units = 8
        left = [1,5,1,4]
        bottom = [5,4,1,1]
        width = [3,2,2,3]
        height = [2,3,3,2]
        
        figsize, coords = PT.make_coordinates(inch_per_unit, x_units, y_units, left, bottom, width, height)

        r = 0.125
        c = [(r, 5*r, 3*r, 2*r), (5*r, 4*r, 2*r, 3*r), (r, r, 2*r, 3*r), (4*r, r, 3*r, 2*r)]
                
        self.assertEqual(figsize, (8,8))
        self.assertEqual(coords, c)


    def test_three_plots(self):
        """
        Three plots, but width has only two elements. 
        
        01234567
        x x  x x

        """
        inch_per_unit = 1.0
        x_units = 8
        y_units = 1
        left = [1,3,6]
        bottom = 0
        width = [1,2]
        height = 1
        
        figsize, coords = PT.make_coordinates(inch_per_unit, x_units, y_units, left, bottom, width, height)
    
        r = 0.125
        c = [(r, 0, r, 1), (3*r, 0, 2*r, 1), (6*r, 0, r, 1)]
                
        self.assertEqual(figsize, (8,1))
        self.assertEqual(coords, c)













 
            
 
class Test_find_longest_list(unittest.TestCase):
    """
    Find the longest list for a series of lists
    
    CHANGELOG:
    20130317/RB: started the suite

    """

    def setUp(self):
        """
        Toggle verbose using the command line. 

        """
        self.flag_verbose = args.verbose

        self.a = [1]
        self.b = [1,2]
        self.c = [1,2,3]
        
        self.d = numpy.array([1,2])
        self.e = numpy.array([1,2,3,4])

    def test_1(self):
        result = PT.find_longest_list(self.a, self.b, self.c)
        self.assertEqual(result, 3) 

    def test_2(self):
        result = PT.find_longest_list(self.a, self.c, self.d)
        self.assertEqual(result, 3)           

    def test_3(self):
        result = PT.find_longest_list(self.a, self.b, self.e)
        self.assertEqual(result, 4)     

    def test_4(self):
        result = PT.find_longest_list(self.a)
        self.assertEqual(result, 1)   
    
if __name__ == '__main__': 

    execute(args)    

























