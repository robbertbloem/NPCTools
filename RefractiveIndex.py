"""
The intention of this script is to process data from RefractiveIndex.info. The data comes in various forms. 
Some data is tabulated. The refractive index values are interpolated for the wavelengths of your choosing. It is not possible to calculate the GVD. Other data comes as a coefficients for equations. The best known is the Sellmeier equation, but there are others (see https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf for the list). This script implements the first six. (I could implement the others, but I didn't find examples of them). 

WARNING: RefractiveIndex.info is a widely used source, but comes with no warrenty. 
- Check the references that are listed.
- The references (and sometimes the comments) contain details about the measurement circumstances, for example temperature, pressure, atmosphere and how the sample is made (thin layer). 
- The data is given for a range, the script will not calculate data outside this range. 

The important functions are get_ri() and get_gvd(). You enter the path to a file and it will import the data. 

The functions plot_ri() and plot_gvd() are wrappers around the get_xx() functions and make a plot. 

The functions ri_for_wavelengths() and gvd_for_wavelengths() are helper functions and assume data in a particular format. 

"""




from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt



import NPCTools.Resources.RI_read_yaml as RIRY
import NPCTools.Resources.Constants as CONST
import NPCTools.Resources.Equations as EQ
import NPCTools.Resources.Mathematics as MATH
import NPCTools.Resources.CommonFunctions as CF
import NPCTools.Debug as DEBUG





def ri_for_wavelengths(db_record, wl_um, interpolate_kind = "default", verbose = 0):
    """
    Calculate the refractive index for wavelengths wl_um. The input is the data from refractiveindex.info and comes as one of 9 equations (not all of them are implemented) or as a table of values. For the latter interpolation will be used to get the values for the asked wavelengths. 
    
    db_record is the dictionary.
      
    """

    if "type" not in db_record:
        raise KeyError("ri_for_wavelengths(): The database record does not have a key 'type'.")
        
    wl_um = CF.make_numpy_ndarray(wl_um)
   
    error_string = "Calculation of the refractive index is not implemented for "
    
    # check the range 
    if "formula" in db_record["type"]:
        
        if wl_um[0] < db_record["range"][0]: 
            raise ValueError("Error, wavelength %1.2f micron is too low! It should be above %1.2f micron." % (wl_um[0], db_record["range"][0]) )
        elif wl_um[-1] > db_record["range"][1]:
            raise ValueError("Error, wavelength %1.2f micron is too high! It should be below %1.2f micron." % (wl_um[-1], db_record["range"][-1]) )
  
    if "tabulated" in db_record["type"]:
        if wl_um[0] < db_record["data"][:,0][0]: 
            raise ValueError("Error, wavelength %1.2f micron is too low! It should be above %1.2f micron." % (wl_um[0], db_record["data"][:,0][0]) )
        elif wl_um[-1] > db_record["data"][:,0][-1]:
            raise ValueError("Error, wavelength %1.2f micron is too high! It should be below %1.2f micron." % (wl_um[-1], db_record["data"][:,0][-1]) )   

    if db_record["type"] == "formula 1":
        DEBUG.verbose("  Using formula 1 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.ones(len(wl_um)) + db_record["coefficients"][0] 
        l = wl_um**2
        for i in range(n_terms):
            ri += (db_record["coefficients"][2*i+1] * l) / (l - db_record["coefficients"][2*i+2]**2)
        ri = numpy.sqrt(ri)
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])

    elif db_record["type"] == "formula 2":
        DEBUG.verbose("  Using formula 2 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.ones(len(wl_um)) + db_record["coefficients"][0] 
        l = wl_um**2
        for i in range(n_terms):
            ri += (db_record["coefficients"][2*i+1] * l) / (l - db_record["coefficients"][2*i+2])
        ri = numpy.sqrt(ri)

        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])

    elif db_record["type"] == "formula 3":
        DEBUG.verbose("  Using formula 3 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.zeros(len(wl_um)) + db_record["coefficients"][0] 
        for i in range(n_terms):
            ri += db_record["coefficients"][2*i+1] * wl_um**db_record["coefficients"][2*i+2]
        ri = numpy.sqrt(ri)
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])
    
    elif db_record["type"] == "formula 4":  
        DEBUG.verbose("  Using formula 4 to calculate refractive indices", verbose_level = 1)
        ri = numpy.zeros(len(wl_um)) + db_record["coefficients"][0] 
        n_coeff = len(db_record["coefficients"])
        if n_coeff >= 5:
            ri += (db_record["coefficients"][1] * wl_um**db_record["coefficients"][2]) / ( wl_um**2  - db_record["coefficients"][3]**db_record["coefficients"][4])
        if n_coeff >= 9:
            ri += (db_record["coefficients"][5] * wl_um**db_record["coefficients"][6]) / ( wl_um**2  - db_record["coefficients"][7]**db_record["coefficients"][8])
        n_terms = int((len(db_record["coefficients"]) - 9)/ 2)        
        for i in range(n_terms):
            ri += db_record["coefficients"][2*i+9] * wl_um**db_record["coefficients"][2*i+10]
        ri = numpy.sqrt(ri)
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])
    
    elif db_record["type"] == "formula 5":
        DEBUG.verbose("  Using formula 5 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.zeros(len(wl_um)) + db_record["coefficients"][0] 
        for i in range(n_terms):
            ri += db_record["coefficients"][2*i+1] * wl_um**db_record["coefficients"][2*i+2]
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])
        
    elif db_record["type"] == "formula 6":
        DEBUG.verbose("  Using formula 6 to calculate refractive indices", verbose_level = 1)
        n_terms = int(len(db_record["coefficients"]) / 2 - 0.5)        
        ri = numpy.ones(len(wl_um)) + db_record["coefficients"][0] 
        l = wl_um**-2
        for i in range(n_terms):
            ri += (db_record["coefficients"][2*i+1]) / (db_record["coefficients"][2*i+2] - l)
        if verbose >= 1:
            for i in range(len(wl_um)):
                print(i, wl_um[i], ri[i])
        
    elif db_record["type"] == "formula 7":
        raise NotImplementedError(error_string + "formula 7")  
        
    elif db_record["type"] == "formula 8":
        raise NotImplementedError(error_string + "formula 8")  
        
    elif db_record["type"] == "formula 9":
        raise NotImplementedError(error_string + "formula 9")  

    elif db_record["type"] == "tabulated n":
        DEBUG.verbose("  Tabulated data (type n)", verbose_level = 1)
        ri = MATH.interpolate_data(db_record["data"][:,0], db_record["data"][:,1], wl_um)
        
    elif db_record["type"] == "tabulated nk":
        DEBUG.verbose("  Tabulated data (type nk)", verbose_level = 1)
        ri = MATH.interpolate_data(db_record["data"][:,0], db_record["data"][:,1], wl_um)
        
    else:
        raise ValueError("The type of data in the database record is unknown (usually formula 1-9 or tabulated data). Type here is %s." % db_record["type"]) 
    
    return ri    


def gvd_for_wavelengths(db_record, wl_um, interpolate_kind = "default", verbose = 0):
    """
    Calculate the GVD for wavelengths wl_um. The GVD is calculated as the second derivative of the refractive index with respect to the wavelength. The input is the data from refractiveindex.info and comes as one of 9 equations (not all of them are implemented) or as a table of values. 
    
    For the equations, the second derivatives were calculated using WolframAlpha. 
    """
    
    if "type" not in db_record:
        raise KeyError("gvd_for_wavelengths: The database record does not have a key 'type'.")
        
    wl_um = CF.make_numpy_ndarray(wl_um)
    
    error_string = "GVD is not implemented for "
    
    # check the range
    if "formula" in db_record["type"]:   
        if wl_um[0] < db_record["range"][0]: 
            raise ValueError("Error, wavelength %1.2f micron is too low! It should be above %1.2f micron." % (wl_um[0], db_record["range"][0]) )
        elif wl_um[-1] > db_record["range"][1]:
            raise ValueError("Error, wavelength %1.2f micron is too high! It should be below %1.2f micron." % (wl_um[-1], db_record["range"][-1]) )
    
    if "tabulated" in db_record["type"]:   
        raise NotImplementedError("GVD can't be calculated for tabulated data.")
        
        
    if db_record["type"] == "formula 1":
        DEBUG.verbose("  Using formula 1 to calculate group velocity dispersion", verbose_level = 1)
        gvd = EQ.gvd_formula_1(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)

    elif db_record["type"] == "formula 2":
        """
        Formula 2 is the same as formula 1, but some given coefficients are already squared. The square root of these coefficients is taken and the GVD equation for formula 1 is used. 
        """
        DEBUG.verbose("  Using formula 2 to calculate group velocity dispersion", verbose_level = 1)
        s = numpy.copy(db_record["coefficients"])
        for i in range(len(s)):
            if i > 0 and i % 2 == 0:
                s[i] = numpy.sqrt(s[i])
        gvd = EQ.gvd_formula_1(wl_um, s)
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)

    elif db_record["type"] == "formula 3":
        DEBUG.verbose("  Using formula 3 to calculate group velocity dispersion", verbose_level = 1)
        gvd = EQ.gvd_formula_3(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)
        
    elif db_record["type"] == "formula 4":
        DEBUG.verbose("  Using formula 4 to calculate group velocity dispersion", verbose_level = 1)
        gvd = EQ.gvd_formula_4(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)
        
    elif db_record["type"] == "formula 5":
        DEBUG.verbose("  Using formula 5 to calculate group velocity dispersion", verbose_level = 1)
        gvd = EQ.gvd_formula_5(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)  
        
    elif db_record["type"] == "formula 6":
        DEBUG.verbose("  Using formula 6 to calculate group velocity dispersion", verbose_level = 1)
        gvd = EQ.gvd_formula_6(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)
        
    elif db_record["type"] == "formula 7":
        DEBUG.verbose("  Using formula 7 to calculate group velocity dispersion", verbose_level = 1)
        gvd = EQ.gvd_formula_7(wl_um, db_record["coefficients"])
        gvd = (1e21 * gvd * wl_um**3) / (2 * numpy.pi * (CONST.c_ms)**2)
        
    elif db_record["type"] == "formula 8":
        raise NotImplementedError(error_string + "formula 8")

    elif db_record["type"] == "formula 9":
        raise NotImplementedError(error_string + "formula 9")

    else:
        raise ValueError("The type of data in the database record is unknown (usually formula 1-9). Type here is %s." % db_record["type"])

    return wl_um, gvd
             
        
def get_ri(paf, wl_um = [], um_range = [0.3, 0.6], n_steps = 100, interpolate_kind = "default", verbose = 0):
    """
    Import a file from RefractiveIndex.info and output the refractive index for the desired wavelengths. 
    
    
    INPUT:
    paf (str): path and filename
    wl_um (ndarray): wavelength axis, in micrometer
    um_range (list with 2 elements): if wl_um is not given, it will plot a range.
    n_steps (int): number of steps to plot the range.
    ax (plt axis): if False, it will make a new figure
    interpolate_kind (str): for tabulated data, the type of interpolation
    
    """
    temp = paf.split("/")
    DEBUG.verbose("Importing data for %s by %s" % (temp[-2], temp[-1][:-4]), verbose_level = 1)

    db_record = RIRY.import_refractive_index(paf = paf, verbose = verbose)
    DEBUG.verbose("  Imported data", verbose_level = 0)
    
    if len(wl_um) == 0:
        wl_um = numpy.linspace(um_range[0], um_range[1], num = n_steps)
    ri = ri_for_wavelengths(db_record, wl_um, verbose = verbose)

    if type(ri) == int:
        return 0, 0

    return wl_um, ri




def effect_gvd(wl_um, gvd, t_fs, d_mm):
    """
    Calculates the effect the group velocity dispersion has on an unchirped Gaussian pulse. 
    
    Source:
    http://www.rp-photonics.com/chromatic_dispersion.html
    numpy.log() is natural log
    
    CHANGELOG:
    20170315/RB: started function. 
    
    INPUT:
    wl_um (ndarray): wavelengths in micron 
    gvd: group velocity dispersion in fs^2/mm
    t: pulse durations
    d: thicknesses of material
    
    OUTPUT:
    t_out (ndarray): 3D array with pulse lengths wl_um x t x d
    
    """
    
    t_fs = CF.make_numpy_ndarray(t_fs)
    d_mm = CF.make_numpy_ndarray(d_mm)

    W, T, D = numpy.meshgrid(wl_um, t_fs, d_mm, indexing = "ij")
    G, dump, dump = numpy.meshgrid(gvd, t_fs, d_mm, indexing = "ij")

    t_out = T * numpy.sqrt(1 + (4 * numpy.log(2) * G * D/ T**2)**2 )
    
    return t_out





def get_gvd(paf, wl_um = [], um_range = [0.3, 0.6], n_steps = 100, interpolate_kind = "default", verbose = 0):
    """
    
    INPUT:
    paf (str): path and filename
    wl_um (ndarray): wavelength axis, in micrometer
    um_range (list with 2 elements): if wl_um is not given, it will plot a range.
    n_steps (int): number of steps to plot the range.
    ax (plt axis): if False, it will make a new figure
    interpolate_kind (str): for tabulated data, the type of interpolation
    
    """
    temp = paf.split("/")
    DEBUG.verbose("Importing data for %s by %s" % (temp[-2], temp[-1][:-4]), verbose_level = 1)

    db_record = RIRY.import_refractive_index(paf = paf, verbose = verbose)
    DEBUG.verbose("  Imported data", verbose_level = 0)

    # error checking
#     del db_record["type"]
#     db_record["type"] = "fiets"

    if len(wl_um) == 0:
        wl_um = numpy.linspace(um_range[0], um_range[1], num = n_steps)
    wl_um, gvd = gvd_for_wavelengths(db_record, wl_um, verbose = verbose)

    if type(gvd) == int:
        return 0, 0

    return wl_um, gvd


def plot_ri(paf, wl_um = [], um_range = [0.3, 0.6], n_steps = 100, ax = False, interpolate_kind = "default", verbose = 0):
    """
    
    INPUT:
    paf (str): path and filename
    wl_um (ndarray): wavelength axis, in micrometer
    um_range (list with 2 elements): if wl_um is not given, it will plot a range.
    n_steps (int): number of steps to plot the range.
    ax (plt axis): if False, it will make a new figure
    interpolate_kind (str): for tabulated data, the type of interpolation
    
    """
    
    wl_um, ri = get_ri(paf = paf, wl_um = wl_um, um_range = um_range, n_steps = n_steps, interpolate_kind = interpolate_kind, verbose = verbose)
        
    if ax == False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    ax.plot(wl_um, ri)

    ax.set_xlabel("Wavelength (micron)")
    ax.set_ylabel("Index of refraction n")
   
    return wl_um, ri



def plot_gvd(paf, wl_um = [], um_range = [0.3, 0.6], n_steps = 100, ax = False, interpolate_kind = "default", verbose = 0):
    """
    
    INPUT:
    paf (str): path and filename
    wl_um (ndarray): wavelength axis, in micrometer
    um_range (list with 2 elements): if wl_um is not given, it will plot a range.
    n_steps (int): number of steps to plot the range.
    ax (plt axis): if False, it will make a new figure
    interpolate_kind (str): for tabulated data, the type of interpolation
    
    """
    
    wl_um, gvd = get_gvd(paf, wl_um = wl_um, um_range = um_range, n_steps = n_steps, interpolate_kind = interpolate_kind, verbose = verbose)

    if ax == False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    ax.plot(wl_um, gvd)
    
    ax.set_xlabel("Wavelength (micron)")
    ax.set_ylabel("GVD (fs^2/mm)")
    
    return wl_um, gvd
      



if __name__ == "__main__": 
    plt.close("all")
    
    path = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndexDB/"
    
#     
#     # formula 1
#     # https://refractiveindex.info/?shelf=main&book=ZnSe&page=Connolly
#     paf = path + "database/main/ZnSe/Connolly.yml"
#     wl_um, ri = get_ri(paf, um_range = [4, 8])
#     wl_um, gvd = get_gvd(paf, um_range = [4, 8])
#     t = [30, 75, 150]
#     d = [1,10]
#     t_out = effect_gvd(wl_um, gvd, t, d)
# #     print(numpy.shape(t_out))
# #     print(gvd[50])
#     _w = 50
#     for _t in range(len(t)):
#         for _d in range(len(d)):
#             plt.plot(wl_um, t_out[:, _t, _d])
#             print(wl_um[_w], t[_t], d[_d], t_out[_w,_t,  _d])
    
#     plt.plot(wl_um, t_out[0,:,0])
#     plt.plot(wl_um, t_out[0,:,1])
    
#     print(t_out)
    
#     a = numpy.arange(3)
#     b = numpy.arange(4) * 10
#     c = numpy.arange(5) * 100
#     
#     A, B = numpy.meshgrid(a, b)
#     
#     print(A)
#     print(B)
#     
#     print(A + B)
#     
#     print(numpy.shape(A), numpy.shape(B)) #, numpy.shape(C))
    
    
#     
#     # formula 2
#     # https://refractiveindex.info/?shelf=main&book=CaF2&page=Daimon-20
#     paf = path + "database/main/CaF2/Daimon-20.yml"
#     wl_um = numpy.array([0.3, 0.4, 0.5, 0.6])
# #     plot_ri(paf, wl_um = wl_um)
#     ri = plot_ri(paf)
#     gvd = plot_gvd(paf, um_range = [0.3, 0.6])
#     print(gvd[0])


    # formula 2
    # https://refractiveindex.info/?shelf=main&book=CaF2&page=Daimon-20
#     paf = path + "database/main/CaF2/Daimon-20.yml"
#     wl_um, gvd = get_gvd(paf, um_range = [0.2, 2.3])
#     t = [30, 75, 150]
#     d = [1,10]
#     t_out = effect_gvd(wl_um, gvd, t, d)
#     _w = 50
#     for _t in range(len(t)):
#         for _d in range(len(d)):
#             plt.plot(wl_um, t_out[:, _t, _d])
#             print(wl_um[_w], t[_t], d[_d], t_out[_w,_t,  _d])


#     paf = path + "database/main/CaF2/Daimon-20.yml"
#     wl_um = [0.8]
#     wl_um, gvd = get_gvd(paf, wl_um = wl_um)
#     t = numpy.arange(10,150)
#     d = [5]
#     t_out = effect_gvd(wl_um, gvd, t, d)
#     plt.plot(t, t_out[0,:,0])
# 
#     # formula 2 AND table
#     # https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
#     paf = path + "database/glass/schott/N-BK7.yml"
# #     wl_um = numpy.array([0.3, 0.4, 0.5, 0.6])
#     um_range = [0.4, 2.4]
# #     plot_ri(paf, um_range = um_range) 
#     wl_um, ri = plot_ri(paf, um_range = um_range)
#     wl_um, gvd = plot_gvd(paf, um_range = um_range) 
#     print(gvd[0])
#     t = numpy.arange(10,150)
#     d = [5]
#     t_out = effect_gvd(wl_um, gvd, t, d)
#     plt.plot(t, t_out[0,:,0])



    # formula 2 AND table
    # https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
    paf = path + "database/glass/schott/N-BK7.yml"
#     wl_um = numpy.array([0.3, 0.4, 0.5, 0.6])
    wl_um = [0.8]

    wl_um, ri = get_ri(paf, wl_um = wl_um)
    wl_um, gvd = get_gvd(paf, wl_um = wl_um) 
    print(gvd[0])
    t = numpy.arange(10,150)
    d = [10]
    t_out = effect_gvd(wl_um, gvd, t, d)
    plt.plot(t, t_out[0,:,0])
    plt.plot(t, t)

#     
#     
# #     paf = path + "database/main/CaF2/Daimon-20.yml"
# #     wl_um = numpy.array([0.3, 20])
# #     plot_ri(paf, wl_um = wl_um)
# # #     plot_ri(paf)
#     
#     
#     paf = path + "database/main/Al2O3/Boidin.yml"
# #     wl_um = numpy.array([0.3, 0.4, 0.5, 0.6])
# #     plot_ri(paf, wl_um = wl_um)
#     plot_ri(paf, um_range = [0.4, 17])
# 
# 
#     # tabulated nk
#     paf = path + "database/main/Ag/Babar.yml"
#     plot_ri(paf, um_range = [0.21, 12])
# 
# 
    # formula 3
#     paf = path + "database/organic/C6H6 - benzene/Moutzouris.yml"
#     ri = plot_ri(paf, um_range = [0.5, 1.5])
#     gvd = plot_gvd(paf, um_range = [0.5, 1.5])
#     print(ri[0])
#     print(gvd[0])
# 
#     # formula 4
#     # https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Eimerl-o
#     paf = path + "database/main/BaB2O4/Eimerl-o.yml"
#     ri = plot_ri(paf, um_range = [0.25, 1.0])
#     gvd = plot_gvd(paf, um_range = [0.25, 1.0])
#     print(ri[0])
#     print(gvd[0])
# 
#     # formula 4
#     # https://refractiveindex.info/?shelf=main&book=BiB3O6&page=Umemura-α
#     paf = path + "database/main/BiB3O6/Umemura-alpha.yml"
#     ri = plot_ri(paf, um_range = [0.5, 1.0])
#     gvd = plot_gvd(paf, um_range = [0.5, 1.0])
#     print(ri[0])
#     print(gvd[0])

 
#     # formula 5
#     # https://refractiveindex.info/?shelf=organic&book=octane&page=Kerl-293K
#     paf = path + "database/organic/C8H18 - octane/Kerl-293K.yml"
#     plot_ri(paf, um_range = [0.33, 0.64])

# 
#     # formula 6
#     # https://refractiveindex.info/?shelf=main&book=H2&page=Peck
#     paf = path + "database/main/H2/Peck.yml"
#     plot_ri(paf, um_range = [0.17, 1.65])

#     n = 101
#     data = numpy.zeros((n, 2))
#     
#     data[:,0] = numpy.arange(n) / 15
#     data[:,1] = numpy.sin(data[:,0])
#     plt.plot(data[:,0], data[:,1])
# 
#     cut_data = data[::5,:]
#     db_record = {"data": cut_data}
# 
#     
#     ri = interpolate_data(db_record, data[:,0], interpolate_kind = "cubic")
#     plt.plot(data[:,0], ri)


    
    plt.show()
    