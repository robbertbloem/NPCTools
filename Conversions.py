from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt


def power_area(power, area, power_unit, area_unit, unit_out, reprate = 2000):

    if power_unit == "MW":
        power_unit = "megawatt"

    power_unit = power_unit.lower()
    area_unit = area_unit.lower()
    unit_out = unit_out.lower()

    energy = power / reprate
    
    if power_unit in ["w", "watt"]:
        power_W = power
    elif power_unit in ["mw", "milliwatt", "milli"]:
        power_W = power / 1e3
    elif power_unit in ["muw", "uw", "microwatt", "micro"]:
        power_W = power / 1e6
    elif power_unit in ["nw", "nanowatt", "nano"]:
        power_W = power / 1e9
    elif power_unit in ["pw", "picowatt", "pico"]:
        power_W = power / 1e12
    elif power_unit in ["kilowatt", "kilo"]:
        power_W = power * 1e3
    elif power_unit in ["megawatt", "mega"]:
        power_W = power * 1e6

    if unit_out in ["mj/cm2"]:
        pass
        
    
    






if __name__ == "__main__": 
    
    power = 1
    area = 1
    power_unit = "W"
    area_unit = "mm"
    unit_out = "mj/cm"
    reprate = 2000
    
    power_area(power, area, power_unit, area_unit, unit_out, reprate = 2000)