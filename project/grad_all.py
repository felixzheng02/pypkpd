"""
## Function written to match MATLAB function grad_all()

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from matpy.matrix import Matrix
from project.lower_tri import lower_tri


def grad_all(func, select_par, nRow, *argv, subset=None, currentOcc=None, noPopED=False, offdiag=False):
    
    arg_list = argv
    def0 = arg_list[select_par]
    
    if currentOcc is None:
        idx = np.arange(1, len(def0)+1)
        idx0 = np.arange(1, len(def0)+1)
        if subset is not None:
            idx0 = np.cumsum(subset) * subset
            idx = np.arange(1, np.max(idx0)+1)
    
        if type(def0) is Matrix and select_par > 6 and offdiag is False:
            idx0 = np.diag(idx0)
    
        if type(def0) == Matrix and select_par > 6 and offdiag is True:
            tmp = 0 * def0.get_data()
            tmp[lower_tri(tmp)] = idx0
            idx0 = tmp + np.transpose(tmp)
    else:
        idx  = np.arange(1, def0.get_shape(0)+1)
        idx0 = np.zeros(def0.get_shape())
        idx0[:, currentOcc] = idx
  

    pypkpd_db = arg_list[len(arg_list)]
    hlf  = pypkpd_db["settings"]["hlf"]
    grad_all_switch = pypkpd_db["settings"]["grad_all_switch"][0]

    if noPopED == True: 
      arg_list = arg_list[-len(arg_list)]
  
    gradX = zeros(nRow, idx.size)

    #Central approximation
    if grad_all_switch == 1:
        for i in idx:
            arg_list[select_par] = def0 + (idx0 == i) * hlf
            def_plus = eval(str(func) + str(arg_list))
            arg_list[select_par] = def0 - (idx0 == i) * hlf
            def_minus = eval(str(func) + str(arg_list))
            if noPopED == False:
                def_plus = def_plus[0]        
                def_minus = def_minus[0]        
            gradX[:,i] = (def_plus - def_minus) / (2.0 * hlf)
    
    else:
    #Complex approximation
        if grad_all_switch == 0:
            for i in idx:
                arg_list[select_par] = def0 + (idx0 == i) * complex(real = 0, imaginary = hlf)
                def_plus = eval(str(func) + str(arg_list))[0]
##没写！！！！pixel image
                gradX[:,i] = Im(def_plus) / hlf
        else:
            raise Exception("Unknown derivative option for grad_all")

    return gradX
   
