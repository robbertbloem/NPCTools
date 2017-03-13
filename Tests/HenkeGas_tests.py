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
import NPCTools.HenkeGas as HG
import NPCTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

global suite_list
suite_list = [
    "Suite 1: Import gas data",
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")
parser.add_argument("-s1", "--skip1", action = "store_true", help = suite_list[0])


# process
args = parser.parse_args()

# reload
if args.reload:
    reload(HG)
    reload(DEBUG)

def execute(args):
    
    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase( Test_import_gas_data)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping :" + suite_list[0], True)
        


class Test_import_gas_data(unittest.TestCase):

    def setUp(self):
        self.flag_verbose = args.verbose

        self.ar_ev_check = numpy.array([10.0, 10.9, 11.8, 12.7, 13.6, 14.5, 15.4, 16.3, 17.2, 18.1, 19.0, 19.9, 20.8, 21.7, 22.6, 23.5, 24.4, 25.3, 26.2, 27.1, 28.0, 28.9, 29.8, 30.7, 31.6, 32.5, 33.4, 34.3, 35.2, 36.1, 37.0, 37.9, 38.8, 39.7, 40.6, 41.5, 42.4, 43.3, 44.2, 45.1, 46.0, 46.9, 47.8, 48.7, 49.6, 50.5, 51.4, 52.3, 53.2, 54.1, 55.0, 55.9, 56.8, 57.7, 58.6, 59.5, 60.4, 61.3, 62.2, 63.1, 64.0, 64.9, 65.8, 66.7, 67.6, 68.5, 69.4, 70.3, 71.2, 72.1, 73.0, 73.9, 74.8, 75.7, 76.6, 77.5, 78.4, 79.3, 80.2, 81.1, 82.0, 82.9, 83.8, 84.7, 85.6, 86.5, 87.4, 88.3, 89.2, 90.1, 91.0, 91.9, 92.8, 93.7, 94.6, 95.5, 96.4, 97.3, 98.2, 99.1, 100.0])
        
        self.ar_data_check = numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.355167671174, -0.375521948494, -0.394048842435, -0.401045433115, -0.398929154199, -0.396552485185, -0.394177987844, -0.3889355352, -0.381868191938, -0.375192117615, -0.368160675551, -0.355472788048, -0.343298231885, -0.330162650039, -0.312176852658, -0.295583494889, -0.273003986413, -0.248774801159, -0.227282326849, -0.206845463155, -0.186118911582, -0.167925599068, -0.151546074925, -0.123545441038, -0.100978809437, -0.0752298250934, -0.0550783201642, -0.0406055886489, -0.0301416309651, -0.022569070643, -0.0176971008089, -0.0141201562015, -0.012137185347, -0.0107816303993, -0.0102743861335, -0.00993210869544, -0.00967002658096, -0.00947911717105, -0.0093016023335, -0.00941254060379, -0.00954126450998, -0.00966558587398, -0.00981215318447, -0.0100476527328, -0.010318857547, -0.0105902318268, -0.0108662287393, -0.0111602247689, -0.0115079314058, -0.0118603798082, -0.0122175813021, -0.0125795473883, -0.012950764123, -0.0133222984296, -0.0136358878736, -0.0138600190467, -0.0140394073048, -0.0142143822274, -0.0143939169476, -0.0145690347874, -0.0147442232673, -0.0149194824444, -0.0150813229447, -0.0151397801723, -0.0151937476749, -0.0152477218846, -0.0152972041369, -0.0153511912065, -0.0154006852498, -0.0154546851855, -0.0154996902614, -0.0154906888731, -0.0154501849344, -0.015387186315, -0.0153286957785, -0.0152702131185, -0.0152117383328, -0.0151532714193, -0.0150948123758, -0.0150408571652, -0.014986908657, -0.0149329668495, -0.014879031741, -0.01482510333, -0.0147711816146, -0.0147217592561, -0.0146723425212, -0.014618439768, -0.0145690347874, -0.014510654328, -0.0144388122267, -0.0143355600223, -0.0142098948102, -0.0140932382344])  
        
        
    def test_import_existing_gas(self):
        """
        Import an existing gas and compare the data with results that were calculated earlier.
        """
        
        ev, data = HG.import_absorption_data("argon", verbose = self.flag_verbose) 

        self.assertTrue(numpy.all(ev == self.ar_ev_check))
        self.assertTrue(numpy.allclose(data, self.ar_data_check))

    def test_import_existing_gas_capitalisation(self):
        """
        Import an existing gas and compare the data with results that were calculated earlier.
        """
        
        ev, data = HG.import_absorption_data("aRgOn", verbose = self.flag_verbose) 

        self.assertTrue(numpy.all(ev == self.ar_ev_check))
        self.assertTrue(numpy.allclose(data, self.ar_data_check))

   
    def test_import_nonexisting_gas(self):
        """
        Try to import a gas for which there is no entry in the list. Check that the output is [0], [0].
        """
    
        ev, data = HG.import_absorption_data("kryptonite", verbose = self.flag_verbose) 
        ev_check = [0]
        data_check = [0]
        self.assertTrue(numpy.all(ev == ev_check))
        self.assertTrue(numpy.allclose(data, data_check))


    def test_import_nonexisting_gas_check(self):
        """
        Try to import a gas for which there is no entry in the list. Double check that it catches that ev_check is [1]. 
        """    
        ev, data = HG.import_absorption_data("kryptonite", verbose = self.flag_verbose) 
        ev_check = [1]
        data_check = [1]
        self.assertFalse(numpy.all(ev == ev_check))
        self.assertFalse(numpy.allclose(data, data_check))


if __name__ == '__main__': 
    
    execute(args)      
















