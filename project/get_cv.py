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


import re
import warnings
import inspect
import numpy as np
from project.size import size
from project.zeros import zeros
from matpy.matrix import Matrix
from project.blockfinal import get_parnam
from project.diag_matlab import diag_matlab
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


def get_rse(fim, poped_db,*argv):
    bpop = poped_db["parameters"]["bpop"][:,1]
    #bpop=poped_db["parameters"]bpop[,2,drop=F],
    d = poped_db["parameters"]["d"][:,1]
    # d = poped_db["parameters"]d[,2,drop=F],
    docc = poped_db["parameters"]["docc"]
    sigma = poped_db["parameters"]["sigma"]
    use_percent = True
    fim_calc_type = poped_db["settings"]["iFIMCalculationType"]
    prior_fim = poped_db["settings"]["prior_fim"]
    #pseudo_on_fail = False,
    ## update poped_db with options supplied in function
    called_args = locals()
    _, _, _, default_args = inspect.getargspec()
    for i in called_args.keys()[-1]:
        if len(re.search('^poped\\_db\\[\\"[\s]*\\"\\]', str(default_args[i]))) == 1:
        #eval(parse(text=paste(capture.output(default_args[[i]]),"=",called_args[[i]])))
        # if (i %in% c('bpop','d')) {
        #   if (eval(parse(text=paste("dim(",i,")[2]>1"))))
        #     (eval(parse(text=paste(i, "=",i,"[,2]"))))
        # }
            eval(str(default_args[[i]]) + "=" + str(i))

    ## if prior is given in poped_db then add it to the given fim
    if len(prior_fim) != 0 and all(size(prior_fim) == size(fim)):
        fim = fim + prior_fim

    try:
        inv_fim = np.linalg.inv(fim)
    except:
        return None
   

    if inv_fim is None:
        mess = "\n  Could not invert the FIM." + "\n  Is the design adequate to estimate all parameters?"
        eig = np.linalg.eig(fim)["values"][0]
        eig.keys() = get_parnam(poped_db)
        neg_vals:np.ndarray = eig[eig < 0]
        num_neg = neg_vals.size
        if num_neg > 0:
            mess = mess + "\n  Potentially problematic parameters and associated eigenvalues:"
            for i in range(0, num_neg):
                mess = mess + ("\n %12s  %8.7e" %(neg_vals[i].keys(),neg_vals[i]))
        #warning(simpleWarning(mess,call="get_rse()"))
        warnings.warn(mess)
        return (np.repeat(np.nan, len(get_parnam(poped_db))))

    param_vars = diag_matlab(inv_fim)
    returnArgs =  get_cv(param_vars,poped_db) 
    params = returnArgs[1]
    params_rse = returnArgs[2]
    parnam = get_parnam(poped_db)
    ret = params_rse[:,:]
    if use_percent: 
        ret[params!=0] = ret[params!=0]*100
    ret.keys() = parnam
    if any(ret==0):
        zero_ret = ret[ret==0].keys()
        mess = "  The following parameters are not estimable:\n  " + ", ".joint(zero_ret) + "\n  Is the design adequate to estimate all parameters?"
        warnings.warn(mess)
        ##warnings.warn(mess, call. = False)
        ret[ret==0] = np.nan
    
    return ret

