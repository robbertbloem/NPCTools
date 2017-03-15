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

# import NPCTools.HenkeGas as HG
import NPCTools.Resources.Equations as EQ
import NPCTools.Debug as DEBUG


# reload(HG)
reload(EQ)
reload(DEBUG)


class Test_quadratic(unittest.TestCase):

    def setUp(self):
        self.verbose = 0
        
    def test_correct(self):
        A = [0,1]
        t = numpy.arange(5)
        res = EQ.quadratic(A, t)
        check = numpy.array([0,1,2,3,4])
        self.assertTrue(numpy.all(res == check))
        
    def test_A_1_element(self):
        A = [0]
        t = numpy.arange(5)
        res = EQ.quadratic(A, t)
        check = numpy.array([0,0,0,0,0])
        self.assertTrue(numpy.all(res == check))

    def test_A_as_int(self):
        A = 0
        t = numpy.arange(5)
        res = EQ.quadratic(A, t)
        check = numpy.array([0,0,0,0,0])
        self.assertTrue(numpy.all(res == check))

    def test_t_as_int(self):
        A = [0,1]
        t = 0 
        res = EQ.quadratic(A, t)
        check = numpy.array([0,0,0,0,0])
        self.assertTrue(numpy.all(res == check))

    def test_t_as_list(self):
        A = [0,1]
        t = [0,0,0,0,0]
        res = EQ.quadratic(A, t)
        check = numpy.array([0,0,0,0,0])
        self.assertTrue(numpy.all(res == check))



class Test_rb_cos(unittest.TestCase):

    def setUp(self):
        self.verbose = 0
        
    def test_correct(self):
        A = [0,1,1,0]
        t = numpy.arange(5) 
        res = EQ.rb_cos(A, t)
        check = numpy.ones(5)
        self.assertTrue(numpy.all(res == check))



if __name__ == '__main__': 
    
    
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_quadratic)
    unittest.TextTestRunner(verbosity=1).run(suite) 
 
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_rb_cos)
    unittest.TextTestRunner(verbosity=1).run(suite) 



