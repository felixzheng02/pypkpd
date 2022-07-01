"""
## Compute the expected parameter relative standard errors 
## 
## This function  computes the expected relative standard errors of a model given a design and a previously computed
## FIM.
## 
## @param fim A Fisher Information Matrix (FIM).
## @param bpop A vector containing the values of the fixed effects used to compute the \code{fim}. 
## @param d A vector containing the values of the diagonals of the between subject variability Matrix.
## @param use_percent Should RSE be reported as percent? 
## @param prior_fim A prior FIM to be added to the \code{fim}. Should be the same size as the \code{fim}.
## @param ... Additional arguments passed to \code{\link{inv}}. 
## @inheritParams evaluate.fim
## @inheritParams Doptim
## @inheritParams create_poped_database
## 
## @return A named list of RSE values.  If the estimated parameter is assumed to be zero then for that 
## parameter the standard error is returned.
## 
## @family evaluate_design
## 
## @example inst/examples_fcn_doc/examples_evaluate.fim.R
## @example tests/testthat/examples_fcn_doc/examples_evaluate.fim.R
## @export

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from project.zeros import zeros
from matpy.matrix import Matrix
from project.get_all_params import get_all_params
from project.get_unfixed_params import get_unfixed_params

def get_cv(param_vars:np.ndarray, poped_db):
    #Return the RSE,CV of parameters
    ## Author: Andrew Hooker
    params_all =  get_all_params(poped_db)[7]
    
    returnArgs =  get_unfixed_params(poped_db, params_all) 
    params:np.ndarray = returnArgs[7]
    var_derivative = returnArgs[8]
    
    if param_vars.size != params.size:
        raise Exception("Number of unfixed parameters not the same as the size of the FIM,\nno RSE can be computed!\n")
    
    params_cv = zeros(size(param_vars))
    for i in range(0, params.size):
        if params[i] != 0:
            if var_derivative[i] == 1:
                params_cv[i] = np.sqrt(param_vars[i])/abs(params[i])
            else:   #Derivative w$r.t to SD instead of var
                params_cv[i] = np.sqrt(param_vars[i])/np.sqrt(params[i])
        
        else:
            params_cv[i] = np.sqrt(param_vars[i])
    
    return {"params": params, "params_cv": params_cv}




