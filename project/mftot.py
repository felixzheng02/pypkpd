"""
#' Evaluate the Fisher Information Matrix (FIM)
#' 
#' Compute the FIM given specific model(s), parameters, design and methods. 
#' 
#' @param poped_db A PopED database.
#' @param bpop The fixed effects parameter values.  Supplied as a vector.
#' @param d A between subject variability matrix (OMEGA in NONMEM).
#' @param docc A between occasion variability matrix.
#' @param sigma A residual unexplained variability matrix (SIGMA in NONMEM).
#' @param model_switch A matrix that is the same size as xt, specifying which model each sample belongs to.
#' @param ni A vector of the number of samples in each group.
#' @param xt A matrix of sample times.  Each row is a vector of sample times for a group.
#' @param x A matrix for the discrete design variables.  Each row is a group.
#' @param a A matrix of covariates.  Each row is a group.
#' @param groupsize A vector of the number of individuals in each group.
#' 
#' @return As a list:
#' \item{ret}{The FIM}
#' \item{poped_db}{A PopED database}
#' 
#' @seealso For an easier function to use, please see \code{\link{evaluate.fim}}.  
#' @family FIM
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_mftot.R
#' @export
#' @keywords internal
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work


## Author: Caiya Zhang, Yuchen Zheng
"""

import numpy as np
from project.size import size
from project.zeros import zeros
from project.mf_all import mf_all
from project.mf_all_loq import mf_all_loq

def mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped_db,*args):
    m = size(ni)[0]
    s = 0
    for i in range(0,m):
        if ni[i]!=0 and groupsize[i]!=0:
            if len(x) != 0:
                x_i = np.transpose(x[i,:])      
            else:
                x_i =  zeros(0,1)
            
            if len(a) != 0:
                a_i = np.transpose(a[i,:])
            else:
                a_i =  zeros(0,1)
            
        # mf_all = function(model_switch,xt,x,a,bpop,d,sigma,docc,poped_db){
        extra_args = [args]
        if extra_args["loq"] is None and extra_args["uloq"] is None: # no loq
            returnArgs = mf_all(np.transpose(model_switch[i,1:ni[i]]),
                                np.transpose(xt[i,1:ni[i]]),
                                x_i,a_i,bpop,d,sigma,docc,poped_db) 
        else:  # handle LOQ 
            returnArgs = mf_all_loq(np.transpose(model_switch[i,1:ni[i]]),
                                np.transposet(xt[i,1:ni[i]]),
                                x_i,a_i,bpop,d,sigma,docc,poped_db,...) 
        
        if returnArgs is None:
            raise Exception("Unknown FIM-calculation type")
        mf_tmp = returnArgs[[0]]
        poped_db = returnArgs[[1]]
        s = s+groupsize[i]*mf_tmp

    return [s, poped_db]
