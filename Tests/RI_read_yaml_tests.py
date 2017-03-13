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
import NPCTools.Resources.RI_read_yaml as RIRY
import NPCTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

global suite_list
suite_list = [
    "Suite 1: coefficients file",
    "Suite 2: data n file",
    "Suite 3: data nk file",
    "Suite 4: string to n col ndarray",
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")
parser.add_argument("-s1", "--skip1", action = "store_true", help = suite_list[0])
parser.add_argument("-s2", "--skip2", action = "store_true", help = suite_list[1])
parser.add_argument("-s3", "--skip3", action = "store_true", help = suite_list[2])
parser.add_argument("-s4", "--skip4", action = "store_true", help = suite_list[3])


# process
args = parser.parse_args()

# reload
if args.reload:
    reload(RIRY)
    reload(DEBUG)

def execute(args):

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_coefficient_file)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping :" + suite_list[0], True)
        
    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_data_n_file)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping :" + suite_list[1], True)

    if args.skip3 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_data_nk_file)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping :" + suite_list[2], True)
   
    if args.skip4 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_string_to_ncol)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping :" + suite_list[3], True)
# 
#     if args.skip5 == False:
#         suite = unittest.TestLoader().loadTestsFromTestCase(Test_2col)
#         unittest.TextTestRunner(verbosity=1).run(suite)    
#     else:
#         DEBUG.verbose("Skipping :" + suite_list[4], True)
# 
#     if args.skip6 == False:
#         suite = unittest.TestLoader().loadTestsFromTestCase(Test_3col)
#         unittest.TextTestRunner(verbosity=1).run(suite)    
#     else:
#         DEBUG.verbose("Skipping :" + suite_list[5], True)



class Test_coefficient_file(unittest.TestCase):

    def setUp(self):
        self.flag_verbose = args.verbose
        
    def test_import_file(self):
        paf = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndex/database/main/CaF2/Daimon-20.yml"
        out = RIRY.import_refractive_index(paf, verbose = 0)
        
#         for i in out:
#             print(i, "->", out[i])

        check = [0, 0.443749998, 0.00178027854, 0.444930066, 0.00788536061, 0.150133991, 0.0124119491, 8.85319946, 2752.28175]
        self.assertTrue(numpy.allclose(out["coefficients"], check))
        
        check = "formula 2"
        self.assertTrue(out["type"] == check)

        self.assertTrue(out["range"][0] == 0.138)
        self.assertTrue(out["range"][1] == 2.326)
        

class Test_data_n_file(unittest.TestCase):

    def setUp(self):
        self.flag_verbose = args.verbose

    def test_import_file(self):
        paf = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndex/database/main/Al2O3/Boidin.yml"
        out = RIRY.import_refractive_index(paf, verbose = 0)
        
#         for i in out:
#             print(i, "->", out[i])      

        check = [0.3, 1.73756]
        self.assertTrue(numpy.allclose(out["data"][0], check))

        check = "tabulated n"
        self.assertTrue(out["type"] == check)
        
        
class Test_data_nk_file(unittest.TestCase):

    def setUp(self):
        self.flag_verbose = args.verbose

    def test_import_file(self):
        paf = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndex/database/main/Ag/Babar.yml"
        out = RIRY.import_refractive_index(paf, verbose = 0)
        
        check = [0.2066, 1.079, 1.247]
        self.assertTrue(numpy.allclose(out["data"][0], check))

        check = "tabulated nk"
        self.assertTrue(out["type"] == check)


    
    def test_import_file_wrong_type_data(self):
        paf = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndex/database/main/Ag/Babar.yml"
        out = RIRY.import_refractive_index(paf, verbose = 0)
        
        check = [0.21, 1.079, 1.247]
        self.assertFalse(numpy.allclose(out["data"][0], check))

        check = "tabulated n"
        self.assertFalse(out["type"] == check)       



class Test_string_to_ncol(unittest.TestCase):
    """
    The conversion from string to an ndarray is done in CommonFunctions and is tested there. These tests are to make sure that conversion from 1D-array to multidimensional array is done properly. 
    
    """


    def setUp(self):
        self.flag_verbose = args.verbose


    def test_1col(self):
        """
        Most basic use of function.
        """
        data = "0 0.443749998 0.00178027854 0.444930066 0.00788536061 0.150133991 0.0124119491 8.85319946 2752.28175"
        out = RIRY.string_to_cols(data, ncols = 1)  
        check = [0, 0.443749998, 0.00178027854, 0.444930066, 0.00788536061, 0.150133991, 0.0124119491, 8.85319946, 2752.28175]
        self.assertTrue(numpy.allclose(out, check))


    def test_2col(self):
        """
        Most basic use of function.
        """
        data = "0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0"
        out = RIRY.string_to_cols(data, ncols = 2)  
        check = [
            [0.0, 1.0],
            [2.0, 3.0],
            [4.0, 5.0],
            [6.0, 7.0]
        ]
        self.assertTrue(numpy.allclose(out, check))

    def test_3col(self):
        """
        Most basic use of function.
        """
        data = "0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8"
        out = RIRY.string_to_cols(data, ncols = 3)  
        check = [
            [0.0, 1.0, 2.0],
            [3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0]
        ]
        self.assertTrue(numpy.allclose(out, check))





if __name__ == '__main__': 
    
    execute(args)      
















