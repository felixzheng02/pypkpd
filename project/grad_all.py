"""
Function translated automatically using 'matlab.to.r()'

Author: Caiya Zhang, Yuchen Zheng
"""
##!!!!f

import numpy as np
from project.zeros import zeros
from project.size import size
from project.feval import do_call

def grad_all (func, select_par, nRow, *args, subset=None, currentOcc=None, noPopED=False, offdiag=False):
    
    for arg in args:
        for i in range(len(args)):
            arg_list = []
            arg_list[i] = arg
    def0 = arg_list[[select_par]]
    
    if currentOcc is None:

        idx = idx0 = len(def0)
        if subset != None:
            idx0 = np.array([])

            for k in range(0, (np.cumsum(subset)*subset).size):
                for i in range(0, (np.cumsum(subset)*subset).shape[1]): #col
                    for j in range(0, (np.cumsum(subset)*subset).shape[0]): #row
                        idx0[k] = np.cumsum(subset)*subset[j,i]
                        k = k + 1
                        j = j + 1
                    j = 0
                    i = i + 1
            
            idx = np.arange(max(idx0))
    
        if type(def0) == np.ndarray and select_par > 6 and offdiag == False:
            idx0 = np.diag(idx0)
    
        if type(def0) == np.ndarray and select_par > 6 and offdiag == True:
            tmp = 0 * def0
###potential problem: lower.tri to np.tril
            tmp[lower.tri(tmp)] = idx0
            idx0 = tmp + np.transpose(tmp)
    else:
        idx  = np.arange(size(def0,1))
        idx0 = zeros(size(def0))
        idx0[:, currentOcc] = idx
  

    poped_db = arg_list[[len(arg_list)]]
    hlf  = poped_db["settings"]["hlf"]
    grad_all_switch = poped_db["settings"]["grad_all_switch"][1]

    if noPopED == True: 
      arg_list = arg_list[-len(arg_list)]
  
    gradX = zeros(nRow, len(idx))

    #Central approximation
    if grad_all_switch == 1:
        for i in idx:
            arg_list[[select_par]] = def0 + (idx0 == i) * hlf
            def_plus = do_call(func, arg_list)
            arg_list[[select_par]] = def0 - (idx0 == i) * hlf
            def_minus = do_call(func, arg_list)
            if noPopED == False:
                def_plus = def_plus[[1]]        
                def_minus = def_minus[[1]]        
            gradX[:,i] = (def_plus - def_minus) / (2.0 * hlf)
    
    else:
    #Complex approximation
        if grad_all_switch == 0:
            for i in idx:
                arg_list[[select_par]] = def0 + (idx0 == i) * complex(real = 0, imaginary = hlf)
                def_plus = do_call(func, arg_list)[[1]]
##没写！！！！pixel image
                gradX[:,i] = Im(def_plus) / hlf
        else:
            raise Exception("Unknown derivative option for grad_all")

    return gradX
   
