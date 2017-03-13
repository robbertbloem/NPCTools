from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import yaml

import urllib

import NPCTools.Resources.CommonFunctions as CF

def string_to_cols(string, ncols = 1, verbose = 0):
    """
    Data is a string. Convert this to a ndarray ncols. 
    """
    data = CF.string_with_numbers_to_list(string)
    if ncols > 1:
        n = int(len(data) / ncols)
        data = numpy.reshape(data, (n, ncols))
#     elif ncols <1:
#         raise ValueError("")
    return data   
  
  
def import_refractive_index(paf, verbose = 0):

    stream = open(paf, "r")
    docs = yaml.load_all(stream)

    output = {}

    for doc in docs:
        
        if verbose >= 1:
            for k,v in doc.items():
                print(k, "->", v)

        if "coefficients" in doc["DATA"][0]:
            output["coefficients"] = string_to_cols( doc["DATA"][0]["coefficients"], ncols = 1, verbose = verbose)
        
        if "type" in doc["DATA"][0]:
            output["type"] = doc["DATA"][0]["type"]
        
        if "data" in doc["DATA"][0]:       
            if output["type"] == "tabulated nk":
                output["data"] = string_to_cols(doc["DATA"][0]["data"], ncols = 3, verbose = verbose)
            elif output["type"] == "tabulated n":            
                output["data"] = string_to_cols(doc["DATA"][0]["data"], ncols = 2, verbose = verbose)
            else:
                print("Unknown data type")
                
        if "range" in doc["DATA"][0]:
            output["range"] = string_to_cols( doc["DATA"][0]["range"], ncols = 1, verbose = verbose)
        
        if "REFERENCES" in doc:
            output["references"] = doc["REFERENCES"]

        if "COMMENTS" in doc:
            output["comments"] = doc["COMMENTS"]
            
        if "INFO" in doc:
            output["info"] = doc["INFO"]

    return output
    




if __name__ == "__main__": 

    pass
#     
#     # coefficients
#     paf = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndex/database/main/CaF2/Daimon-20.yml"
# #     
#     # tabulated nk
#     paf = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndex/database/main/Ag/Babar.yml"
#     
#     # tabulated n
#     paf = "/Users/rbloem/Developer/NPCTools/Data/RefractiveIndex/database/main/Al2O3/Boidin.yml"
#     
#     
#     import_refractive_index(paf)