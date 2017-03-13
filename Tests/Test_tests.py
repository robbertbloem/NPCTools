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



def execute(args):

#     print(args)

    for i in range(len(suite_list)):

        if args.skip0 == False:
            suite = unittest.TestLoader().loadTestsFromTestCase(Test_suite_1)
            unittest.TextTestRunner(verbosity=1).run(suite)    
        else:
            DEBUG.verbose("Skipping :" + suite_list[i]["name"], True)
        



class Test_suite_1(unittest.TestCase):
    """
    The conversion from string to an ndarray is done in CommonFunctions and is tested there. These tests are to make sure that conversion from 1D-array to multidimensional array is done properly. 
    
    """


    def setUp(self):
        self.flag_verbose = args.verbose


    def test_s1_1(self):
        x = 1 #numpy.array([1,2,3])
        self.assertTrue(x == 1)  #numpy.all(x, numpy.array([1,2,3])))  



class Test_suite_2(unittest.TestCase):
    """
    The conversion from string to an ndarray is done in CommonFunctions and is tested there. These tests are to make sure that conversion from 1D-array to multidimensional array is done properly. 
    
    """


    def setUp(self):
        self.flag_verbose = args.verbose


    def test_s2_1(self):
        x = 1 #numpy.array([1,2,3])
        self.assertTrue(x == 1)  #numpy.all(x, numpy.array([1,2,3])))  



# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

global suite_list
suite_list = [
    {"idx": 0, "name": "Coefficients file", "suite": Test_suite_1, "help": "This is not helping"},
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")

for i in range(len(suite_list)):
    short_name = "-s%i" % i    
    long_name = "--skip%i" % i  
    parser.add_argument(short_name, long_name, action = "store_true", help = suite_list[i]["help"])
    
    
# parser.add_argument("-s2", "--skip2", action = "store_true", help = suite_list[1])
# parser.add_argument("-s3", "--skip3", action = "store_true", help = suite_list[2])
# parser.add_argument("-s4", "--skip4", action = "store_true", help = suite_list[3])


# process
args = parser.parse_args()

# reload
if args.reload:
    reload(RIRY)
    reload(DEBUG)




if __name__ == '__main__': 
    
    execute(args)      
















