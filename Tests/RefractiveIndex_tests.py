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
import NPCTools.RefractiveIndex as RI
import NPCTools.Debug as DEBUG


# reload(HG)
reload(RI)
reload(DEBUG)


class Test_ri_gvd_with_formulas(unittest.TestCase):
    """
    This test checks if the calculated values correspond to the values from refractiveindex.info. It checks this both by importing the files and by calculating it from info from the files. 
    """


    def setUp(self):
        self.verbose = 0

        self.rtol_ri = 0.0001
        self.rtol_gvd = 0.001

        path = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndexDB/"

        # formula 1
        # https://refractiveindex.info/?shelf=main&book=ZnSe&page=Connolly
        self.f1_paf = path + "database/main/ZnSe/Connolly.yml"
        self.f1 = {
            "type": "formula 1", 
            "range": [0.54, 18.2], 
            "coefficients": [0, 4.45813734, 0.200859853, 0.467216334, 0.391371166, 2.89566290, 47.1362108],
        }

        # formula 2
        # https://refractiveindex.info/?shelf=main&book=CaF2&page=Daimon-20
        self.f2_paf = path + "database/main/CaF2/Daimon-20.yml"
        self.f2 = {
            "type": "formula 2", 
            "range": [0.138, 2.326], 
            "coefficients": [0, 0.443749998, 0.00178027854, 0.444930066, 0.00788536061, 0.150133991, 0.0124119491, 8.85319946, 2752.28175],
        }

        # formula 3
        # https://refractiveindex.info/?shelf=organic&book=benzene&page=Moutzouris
        self.f3_paf = path + "database/organic/C6H6 - benzene/Moutzouris.yml"    
        self.f3 = {
            "type": "formula 3", 
            "range": [0.450, 1.551], 
            "coefficients": [2.170184597, 0.00059399, 2, 0.02303464, -2, -0.000499485, -4, 0.000178796, -6],
        }

        # formula 4
        # https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Eimerl-o
        self.f4_paf = path + "database/main/BaB2O4/Eimerl-o.yml"            
        self.f4 = {
            "type": "formula 4", 
            "range": [0.22, 1.06], 
            "coefficients": [2.7405, 0.0184, 0, 0.0179, 1, 0, 0, 0, 1, -0.0155, 2],
        }

        # formula 5
        # https://refractiveindex.info/?shelf=organic&book=octane&page=Kerl-293K
        self.f5_paf = path + "database/organic/C8H18 - octane/Kerl-293K.yml"            
        self.f5 = {
            "type": "formula 5", 
            "range": [0.326, 0.644], 
            "coefficients": [1.39260498, -4.48963e-3, -1, 4.79591e-3, -2],
        }

        # formula 6
        # https://refractiveindex.info/?shelf=main&book=H2&page=Peck
        self.f6_paf = path + "database/main/H2/Peck.yml"           
        self.f6 = {
            "type": "formula 6", 
            "range": [0.1680, 1.6945], 
            "coefficients": [0, 0.0148956, 180.7, 0.0049037, 92],
        }
        

        # tabulated n
        self.tab_n_paf = path + "database/main/Al2O3/Boidin.yml"

        # tabulated nk
        self.tab_nk_paf = path + "database/main/Ag/Babar.yml"
    
    def test_formula_1_ri(self):
        wl_um = numpy.array([0.6, 4.0, 5.0, 6.0, 10.0, 18.0])
        ri = RI.ri_for_wavelengths(self.f1, wl_um, verbose = self.verbose) 
        ri_check = numpy.array([2.6141, 2.4331, 2.4295, 2.4258, 2.4065, 2.3306])
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))
             
        wl_um, ri = RI.get_ri(self.f1_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))


    def test_formula_1_gvd(self):
        wl_um = numpy.array([0.6, 4.0, 5.0, 6.0, 10.0, 18.0])
        gvd = RI.gvd_for_wavelengths(self.f1, wl_um, verbose = self.verbose) 
        gvd_check = numpy.array([2271.1, 75.299, -17.134, -136.21, -1221.0, -14134])
        self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))

        wl_um, gvd = RI.get_gvd(self.f1_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))
       

    def test_formula_2_ri(self):
        wl_um = numpy.array([0.15, 0.4, 0.6, 0.8, 1.5, 2.3])
        ri = RI.ri_for_wavelengths(self.f2, wl_um, verbose = self.verbose) 
        ri_check = numpy.array([1.5817, 1.4419, 1.4336, 1.4306, 1.4263, 1.4223])
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))
        
        wl_um, ri = RI.get_ri(self.f2_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))


    def test_formula_2_gvd(self):
        wl_um = numpy.array([0.15, 0.4, 0.6, 0.8, 1.5, 2.3])
        gvd = RI.gvd_for_wavelengths(self.f2, wl_um, verbose = self.verbose) 
        gvd_check = numpy.array([846.01, 67.513, 40.230, 27.796, 1.8579, -39.704])
        self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))
        
        wl_um, gvd = RI.get_gvd(self.f2_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))


    def test_formula_3_ri(self):
        wl_um = numpy.array([0.5, 1.0, 1.5])
        ri = RI.ri_for_wavelengths(self.f3, wl_um, verbose = self.verbose) 
        ri_check = numpy.array([1.5053, 1.4810, 1.4770])
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))
        
        wl_um, ri = RI.get_ri(self.f3_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))


    def test_formula_3_gvd(self):
        wl_um = numpy.array([0.5, 1.0, 1.5])
        gvd = RI.gvd_for_wavelengths(self.f3, wl_um, verbose = self.verbose) 
        gvd_check = numpy.array([253.85, 81.590, 56.390])
        self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))
        
        wl_um, gvd = RI.get_gvd(self.f3_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))
        

    def test_formula_4_ri(self):
        wl_um = numpy.array([0.25, 0.6, 1.0])
        ri = RI.ri_for_wavelengths(self.f4, wl_um, verbose = self.verbose) 
        ri_check = numpy.array([1.7754, 1.6699, 1.6564])
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))
        
        wl_um, ri = RI.get_ri(self.f4_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))


    def test_formula_4_gvd(self):
        wl_um = numpy.array([0.25, 0.6, 1.0])
        gvd = RI.gvd_for_wavelengths(self.f4, wl_um, verbose = self.verbose) 
        gvd_check = numpy.array([637.15, 111.14, 45.633])
        self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))
        
        wl_um, gvd = RI.get_gvd(self.f4_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))


    def test_formula_5_ri(self):
        wl_um = numpy.array([0.35, 0.5, 0.63])
        ri = RI.ri_for_wavelengths(self.f5, wl_um, verbose = self.verbose) 
        ri_check = numpy.array([1.4189, 1.4028, 1.3976])
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))
        
        wl_um, ri = RI.get_ri(self.f5_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))


    # gvd_check is correct, but this is not implemented yet 
#     def test_formula_5_gvd(self):
#         wl_um = numpy.array([0.35, 0.5, 0.63])
#         gvd = RI.gvd_for_wavelengths(self.f5, wl_um, verbose = self.verbose) 
#         gvd_check = numpy.array([129.69, 86.012, 64.982])
#         self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))
# 
#         wl_um, gvd = RI.get_gvd(self.f5_paf, wl_um = wl_um, verbose = self.verbose)
#         self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))


    def test_formula_6_ri(self):
        wl_um = numpy.array([0.17, 0.7, 1.2, 1.65])
        ri = RI.ri_for_wavelengths(self.f6, wl_um, verbose = self.verbose) 
        ri_check = numpy.array([1.0001874, 1.0001379, 1.0001365, 1.0001361])
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))
        
        wl_um, ri = RI.get_ri(self.f6_paf, wl_um = wl_um, verbose = self.verbose)
        self.assertTrue(numpy.allclose(ri, ri_check, rtol = self.rtol_ri))

    # gvd_check is correct, but this is not implemented yet 
#     def test_formula_6_gvd(self):
#         wl_um = numpy.array([0.17, 0.7, 1.2, 1.65])
#         gvd = RI.gvd_for_wavelengths(self.f6, wl_um, verbose = self.verbose) 
#         gvd_check = numpy.array([0.22527, 0.016515, 0.010617, 0.0081099])
#         self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))
# 
#         wl_um, gvd = RI.get_gvd(self.f6_paf, wl_um = wl_um, verbose = self.verbose)
#         self.assertTrue(numpy.allclose(gvd, gvd_check, rtol = self.rtol_gvd))



class Test_interpolate_data(unittest.TestCase):
    """
    Some basic checks of the interpolation routine. It checks if it works, but doesn't check the data. 
    """


    def setUp(self):
        self.verbose = 0

    def test_default(self):
        
        data = numpy.array([[0.0, 2.0], [1.0, 3.0]])        
        db_record = {"data": data}
        wl_um = numpy.array([0, 0.5, 1.0])
        ri = RI.interpolate_data(db_record, wl_um, interpolate_kind = "default", verbose = self.verbose)
        ri_check = numpy.array([2.0, 2.5, 3.0])
        self.assertTrue(numpy.all(ri == ri_check))
        
    def test_linear_1(self):
        
        data = numpy.array([[0.0, 2.0], [1.0, 3.0]])        
        db_record = {"data": data}
        wl_um = numpy.array([0, 0.5, 1.0])
        ri = RI.interpolate_data(db_record, wl_um, interpolate_kind = "linear", verbose = self.verbose)
        ri_check = numpy.array([2.0, 2.5, 3.0])
        self.assertTrue(numpy.all(ri == ri_check))
        
    def test_linear_2(self):

        n = 101
        data = numpy.zeros((n, 2))   
        data[:,0] = numpy.arange(n) / 15
        data[:,1] = numpy.sin(data[:,0])        
        cut_data = data[::5,:]
        db_record = {"data": cut_data}   
        ri = RI.interpolate_data(db_record, data[:,0], interpolate_kind = "linear", verbose = self.verbose)
        
    def test_cubic_1(self):

        n = 101
        data = numpy.zeros((n, 2))   
        data[:,0] = numpy.arange(n) / 15
        data[:,1] = numpy.sin(data[:,0])        
        cut_data = data[::5,:]
        db_record = {"data": cut_data}   
        ri = RI.interpolate_data(db_record, data[:,0], interpolate_kind = "cubic", verbose = self.verbose)

    def test_quadratic_1(self):

        n = 101
        data = numpy.zeros((n, 2))   
        data[:,0] = numpy.arange(n) / 15
        data[:,1] = numpy.sin(data[:,0])        
        cut_data = data[::5,:]
        db_record = {"data": cut_data}   
        ri = RI.interpolate_data(db_record, data[:,0], interpolate_kind = "quadratic", verbose = self.verbose)




if __name__ == '__main__': 
    
    
    suite = unittest.TestLoader().loadTestsFromTestCase( Test_ri_gvd_with_formulas)
    unittest.TextTestRunner(verbosity=1).run(suite) 
 
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_interpolate_data)
    unittest.TextTestRunner(verbosity=1).run(suite) 
 
 
 



