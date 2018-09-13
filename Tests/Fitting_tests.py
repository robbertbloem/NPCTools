"""

"""

import imp
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import NPCTools.Resources.Fitting as FI

imp.reload(FI)

class Test_basic(unittest.TestCase):

    def setUp(self):
        self.verbose = 0
        
    def test_set_function(self):
        f = FI.Fitting()
        f.set_function(FI.gaussian)
        
        xdata = numpy.arange(-100,100)
        f.set_xdata(xdata)
        
#         mu, sigma, amplitude, offset
        params = (0, 10, 1, 0)
        
        y = f.calculate_from_params(params)
        print(y)

if __name__ == "__main__": 


    suite = unittest.TestLoader().loadTestsFromTestCase(Test_basic)
    unittest.TextTestRunner(verbosity=1).run(suite) 