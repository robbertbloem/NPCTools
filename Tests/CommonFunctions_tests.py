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

# import NPCTools.ClassTools as CT
import NPCTools.Resources.CommonFunctions as CF
import NPCTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

global suite_list
suite_list = [
    "Suite 1: Test make numpy_ndarray",
    "Suite 2: Test string with numbers to list",
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
    reload(CF)
    reload(DEBUG)

def execute(args):
    
    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase( Test_make_numpy_ndarray)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping :" + suite_list[0], True)
        
    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase( Test_string_with_numbers_to_list)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping :" + suite_list[1], True)



   
class Test_make_numpy_ndarray(unittest.TestCase):
    """
    

    CHANGELOG:
    20130317/RB: started the suite

    """

    def setUp(self):
        """
        Toggle verbose using the command line. 
        
        int, float, string: make it a list and then numpy.ndarray.
        list: make it a numpy.ndarray
        tuple: convert to a list, then numpy.ndarray
        numpy.ndarray: return directly
        dict: give error
        
        self.values = [[val, type of return]], if type of return is numpy.ndarray, leave it away.
        
        """

        self.flag_verbose = args.verbose
        
        self.values = [
            [1],                        # int
            [0.2],                      # float
            ["fiets"],                  # string
            [(1,2)],                    # tuple
            [[1]],                      # short list
            [[1,2]],                    # longer list
            [numpy.array([1,2])],       # numpy.ndarray
        ]

    def test_int_1(self):
        result = CF.make_numpy_ndarray(1)
        self.assertEqual(type(result), numpy.ndarray)  
        self.assertEqual(result, numpy.array(1))
        
    def test_float_1(self):
        result = CF.make_numpy_ndarray(1)
        self.assertEqual(type(result), numpy.ndarray)  
        self.assertEqual(result, numpy.array(1))

    def test_string_1(self):
        result = CF.make_numpy_ndarray("fiets")
        self.assertEqual(type(result), numpy.ndarray)  
        self.assertEqual(result, ["fiets"])
        self.assertTrue(result[0] == "fiets")

    def test_string_2(self):
        result = CF.make_numpy_ndarray(["fiets", "auto"])
        self.assertEqual(type(result), numpy.ndarray)  
        self.assertTrue(numpy.all(numpy.array(["fiets", "auto"]) == result)) 
        
    def test_tuple_1(self):
        result = CF.make_numpy_ndarray((1,2))
        self.assertEqual(type(result), numpy.ndarray)  
        self.assertTrue(numpy.all(numpy.array([1,2]) == result)) 

    def test_list_1(self):
        result = CF.make_numpy_ndarray([1,2])
        self.assertEqual(type(result), numpy.ndarray)  
        self.assertTrue(numpy.all(numpy.array([1,2]) == result)) 

    
    @unittest.expectedFailure     
    def test_dict(self):
        result = CF.make_numpy_ndarray({"a": 1, "b": 2})
        


class Test_string_with_numbers_to_list(unittest.TestCase):
    """
    CHANGELOG:
    20170309/RB: started the suite

    """

    def setUp(self):

        self.flag_verbose = args.verbose
    
    
    def test_string_with_float_commas(self):      
        string = "0.0, 1.1, 2.2"
        check = [0.0, 1.1, 2.2]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))


    def test_string_with_float_commas_enter(self):      
        string = "0.0, 1.1\n 2.2"
        check = [0.0, 1.1, 2.2]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))


    def test_string_with_float_spaces(self):        
        string = "0.0 1.1 2.2"
        check = [0.0, 1.1, 2.2]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))


    def test_string_with_float_spaces_enter(self):       
        string = "0.0 1.1\n 2.2"
        check = [0.0, 1.1, 2.2]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))

        
    def test_string_with_int_commas(self):      
        string = "0, 1, 2"
        check = [0.0, 1.0, 2.0]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))


    def test_string_with_int_commas_enter(self):      
        string = "0, 1\n 2"
        check = [0.0, 1.0, 2.0]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))


    def test_string_with_int_spaces(self):        
        string = "0 1 2"
        check = [0.0, 1.0, 2.0]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))


    def test_string_with_int_spaces_enter(self):       
        string = "0 1\n 2"
        check = [0, 1, 2]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))      
        
        
    def test_string_with_int_to_float(self):       
        string = "0 1\n 2"
        check = [0.0, 1.0, 2.0]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))      


    def test_string_with_exp_to_float(self):      
        string = "1e+3, 0.5 1.2e-3"
        check = [1000.0, 0.5, 0.0012]
        res = CF.string_with_numbers_to_list(string)
        self.assertTrue(numpy.allclose(res, check))




if __name__ == '__main__': 
    
    execute(args)      
















