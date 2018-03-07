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


class Test_reflectance(unittest.TestCase):

    def setUp(self):
        self.verbose = 0
        
    def test_a_deg_simple(self):
        """
        Basic test using a_deg
        """
        n = 10
        n1 = numpy.ones(n)
        n2 = numpy.linspace(1.1, 1.5, num = n)
        a_deg = [0, 10]
        a_deg, Rs, Rp = EQ.reflectance(n1, n2, a_deg = a_deg)
        self.assertTrue(numpy.shape(a_deg) == (2,))
        self.assertTrue(numpy.shape(Rs) == (2,10))
        self.assertTrue(numpy.shape(Rp) == (2,10))
        
    def test_arange_simple(self):
        """
        Basic test using a_range
        """
        n = 10
        n1 = numpy.ones(n)
        n2 = numpy.linspace(1.1, 1.5, num = n)
        a_range = (0, 90)
        a_deg, Rs, Rp = EQ.reflectance(n1, n2, a_range = a_range)
        self.assertTrue(numpy.shape(a_deg) == (91,))
        self.assertTrue(numpy.shape(Rs) == (91,10))
        self.assertTrue(numpy.shape(Rp) == (91,10))
        self.assertTrue(numpy.count_nonzero(numpy.isnan(Rs)) == 0)


    def test_arange_critical_angle(self):
        """
        Check if the critical angle is handled properly (it should give NaN)
        """
        n = 10
        n2 = numpy.ones(n)
        n1 = numpy.linspace(1.1, 1.5, num = n)
        a_range = (0, 90)
        a_deg, Rs, Rp = EQ.reflectance(n1, n2, a_range = a_range)
#         print(numpy.count_nonzero(numpy.isnan(Rs)))
        self.assertTrue(numpy.shape(a_deg) == (91,))
        self.assertTrue(numpy.shape(Rs) == (91,10))
        self.assertTrue(numpy.shape(Rp) == (91,10))       
        self.assertTrue(numpy.count_nonzero(numpy.isnan(Rs)) == 388)

    def test_arange_critical_angle_2(self):
        """
        Check if the critical angle is handled properly (it should give NaN)
        """
        n = 10
        n1 = numpy.ones(n)
        n2 = numpy.linspace(0.5, 1.5, num = n)
        a_range = (0, 90)
        a_deg, Rs, Rp = EQ.reflectance(n1, n2, a_range = a_range)
#         print(numpy.count_nonzero(numpy.isnan(Rs)))
        self.assertTrue(numpy.shape(a_deg) == (91,))
        self.assertTrue(numpy.shape(Rs) == (91,10))
        self.assertTrue(numpy.shape(Rp) == (91,10))       
        self.assertTrue(numpy.count_nonzero(numpy.isnan(Rs)) == 211)



    def test_arange_all_critical_angle(self):
        """
        Check if the critical angle is handled properly (it should give NaN)
        This tests for between 50 and 90 degrees, all should be NaN
        """
        n = 10
        n2 = numpy.ones(n)
        n1 = numpy.linspace(1.1, 1.5, num = n)
        a_range = (70, 90)
        a_deg, Rs, Rp = EQ.reflectance(n1, n2, a_range = a_range)
#         print(numpy.count_nonzero(numpy.isnan(Rs)))
        self.assertTrue(numpy.shape(a_deg) == (21,))
        self.assertTrue(numpy.shape(Rs) == (21,10))
        self.assertTrue(numpy.shape(Rp) == (21,10))       
        self.assertTrue(numpy.count_nonzero(numpy.isnan(Rs)) == 210)
#         print(Rs)

    def test_arange_no_critical_angle(self):
        """
        Check if the critical angle is handled properly (it should give NaN)
        This tests for between 0 and 40 degrees, none should be NaN       
        """
        n = 10
        n2 = numpy.ones(n)
        n1 = numpy.linspace(1.1, 1.5, num = n)
        a_range = (0, 40)
        a_deg, Rs, Rp = EQ.reflectance(n1, n2, a_range = a_range)
#         print(numpy.count_nonzero(numpy.isnan(Rs)))
        self.assertTrue(numpy.shape(a_deg) == (41,))
        self.assertTrue(numpy.shape(Rs) == (41,10))
        self.assertTrue(numpy.shape(Rp) == (41,10))       
        self.assertTrue(numpy.count_nonzero(numpy.isnan(Rs)) == 0)


    def test_int_inputs(self):
        """
        use integer inputs for a_deg, n1 and n2
        """
        n1 = 0 
        n2 = 1.5 
        a_deg = 0
        a_deg, Rs, Rp = EQ.reflectance(n1, n2, a_deg = a_deg)
#         print(numpy.count_nonzero(numpy.isnan(Rs)))
        self.assertTrue(numpy.shape(a_deg) == (1,))
        self.assertTrue(numpy.shape(Rs) == (1,1))
        self.assertTrue(numpy.shape(Rp) == (1,1))

    def test_list_inputs(self):
        """
        use integer inputs for a_deg, n1 and n2
        """
        n1 = [1,1,1,1,1]
        n2 = [1.1, 1.2, 1.3, 1.4, 1.5]
        a_deg = [0, 10]
        a_deg, Rs, Rp = EQ.reflectance(n1, n2, a_deg = a_deg)
#         print(numpy.count_nonzero(numpy.isnan(Rs)))
        self.assertTrue(numpy.shape(a_deg) == (2,))
        self.assertTrue(numpy.shape(Rs) == (2,5))
        self.assertTrue(numpy.shape(Rp) == (2,5))


if __name__ == '__main__': 
    
    
#     suite = unittest.TestLoader().loadTestsFromTestCase(Test_quadratic)
#     unittest.TextTestRunner(verbosity=1).run(suite) 
#  
#     suite = unittest.TestLoader().loadTestsFromTestCase(Test_rb_cos)
#     unittest.TextTestRunner(verbosity=1).run(suite) 

    suite = unittest.TestLoader().loadTestsFromTestCase(Test_reflectance)
    unittest.TextTestRunner(verbosity=1).run(suite) 


