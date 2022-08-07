"""
## Evaluate the Fisher Information Matrix (FIM)
## 
## Compute the FIM given specific model(s), parameters, design and methods. 
## 
## @param pypkpd_db A PopED database.
## @param bpop The fixed effects parameter values.  Supplied as a vector.
## @param d A between subject variability Matrix (OMEGA in NONMEM).
## @param docc A between occasion variability Matrix.
## @param sigma A residual unexplained variability Matrix (SIGMA in NONMEM).
## @param model_switch A Matrix that is the same size as xt, specifying which model each sample belongs to.
## @param ni A vector of the number of samples in each group.
## @param xt A Matrix of sample times.  Each row is a vector of sample times for a group.
## @param x A Matrix for the discrete design variables.  Each row is a group.
## @param a A Matrix of covariates.  Each row is a group.
## @param groupsize A vector of the number of individuals in each group.
## 
## @return As a list:
## \item{ret}{The FIM}
## \item{pypkpd_db}{A PopED database}
## 
## @seealso For an easier function to use, please see \code{\link{evaluate.fim}}.  
## @family FIM
## 
## @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
## @example tests/testthat/examples_fcn_doc/examples_mftot.R
## @export
## @keywords internal
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work


## Author: Caiya Zhang, Yuchen Zheng
"""

import numpy as np
from project.size import size
from project.mf_all import mf_all
from project.mf_all_loq import mf_all_loq

def mftot(model_switch,
            groupsize,
            ni,
            xt,
            x,
            a,
            bpop,
            d,
            sigma,
            docc,
            pypkpd_db,
            **kwargs):

    m = ni.get_shape(0)
    s = 0
    for i in range(0, m):
        if ni.get_data()[0, i] != 0 and groupsize.get_data()[0, i] != 0:
            if x.get_size() != 0:
                x_i = np.transpose(x.get_data()[i, :])      
            else:
                x_i =  np.zeros([0, 1])
            
            if a.get_size() != 0:
                a_i = np.transpose(a[i, :])
            else:
                a_i =  np.zeros([0, 1])
            
        # mf_all = function(model_switch,xt,x,a,bpop,d,sigma,docc,pypkpd_db){
        extra_args = kwargs
        if "log" not in extra_args.keys() and "ulog" not in extra_args.keys(): # no loq
            returnArgs = mf_all(np.transpose(model_switch.get_data()[i, 0:ni.get_data()[0, i]]),
                                np.transpose(xt.get_data()[i, 0:ni.get_data()[0, i]]),
                                x_i, a_i, bpop, d, sigma, docc, pypkpd_db) 
        else:  # handle LOQ 
            returnArgs = mf_all_loq(np.transpose(model_switch.get_data()[i, 0:ni.get_data()[0, i]]),
                                np.transpose(xt.get_data()[i, 0:ni.get_data()[0, i]]),
                                x_i, a_i, bpop, d, sigma, docc, pypkpd_db, kwargs) 
        
        if returnArgs is None:
            raise Exception("Unknown FIM-calculation type")
        mf_tmp = returnArgs[0]
        pypkpd_db = returnArgs[1]
        s = s + groupsize.get_data()[0, i] * mf_tmp

    return {"ret": s, "pypkpd_db": pypkpd_db}
