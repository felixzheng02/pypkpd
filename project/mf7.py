"""
#' The full Fisher Information Matrix (FIM) for one individual Calculating one model switch at a time, good for large matrices.
#' 
#' Compute the full FIM for one individual given specific model(s), parameters, design and methods. 
#' This computation calculates the FIM for each model switch separately.  Correlations between the models parameters are assumed to be zero.
#' 
#' @inheritParams mf3
#' 
#' @return As a list:
#' \item{ret}{The FIM for one individual}
#' \item{poped_db}{A PopED database}
#' 
#' @family FIM
#' 
# @example tests/testthat/examples_fcn_doc/warfarin_basic.R
# @example tests/testthat/examples_fcn_doc/examples_mf7.R
# @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'


## Author: Caiya Zhang, Yuchen Zheng
"""

import numpy as np
from project.m1 import m1
from project.m2 import m2
from project.m3 import m3
from project.size import size
from project.zeros import zeros

def mf7(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped_db):
  
    #This calculation of FIM divides the calculation up into one calculation
    #per model switch
    
    numnotfixed_bpop = sum(poped_db["parameters"]["notfixed_bpop"])
    numnotfixed_d    = sum(poped_db["parameters"]["notfixed_d"])
    numnotfixed_covd = sum(poped_db["parameters"]["notfixed_covd"])
    numnotfixed_docc  = sum(poped_db["parameters"]["notfixed_docc"])
    numnotfixed_covdocc  = sum(poped_db["parameters"]["notfixed_covdocc"])
    numnotfixed_sigma  = sum(poped_db["parameters"]["notfixed_sigma"])
    numnotfixed_covsigma  = sum(poped_db["parameters"]["notfixed_covsigma"])
    
    ret = 0
    
    for i in range(0,poped_db["settings"]["iFOCENumInd"]):
        b_ind = poped_db["parameters"]["b_global"][:,i]
        bocc_ind = poped_db["parameters"]["bocc_global"][[i]]
        
        for j in (0,max(model_switch)):
            xt_new = zeros(sum(model_switch==j),1)
            m = 1
            for k in range(0,model_switch.size):
                if model_switch[k] == j:
                    xt_new[m]=xt_ind[k]
                    m=m+1

        n = size(xt_new)[0]
        model_switch_new = np.arrary(1,sum(model_switch==j),1)*j
        
        if len(xt_new) != 0:
            f1 = zeros(n+n*n,numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)
            returnArgs = m1(model_switch_new,xt_new,x,a,bpop,b_ind,bocc_ind,d,sigma,poped_db) 
            f1[1:n,1:numnotfixed_bpop] = returnArgs[1]
            poped_db = returnArgs[[2]]
            returnArgs = m2(model_switch_new,xt_new,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db) 
            f1[(n+1):(n+n*n),1:numnotfixed_bpop] = returnArgs[1]
            poped_db = returnArgs[[2]]
            returnArgs = m3(model_switch_new,xt_new,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,True,poped_db) 
            f1[(n+1):(n+n*n),(numnotfixed_bpop+1):(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)] = returnArgs[0]
            poped_db = returnArgs[[2]]
            f2 = zeros(n+n*n,n+n*n)
            returnArgs =  v(model_switch_new,xt_new,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db) 
            v_tmp = returnArgs[1]
            poped_db = returnArgs[[2]]
            if any(v_tmp != 0):
            
                v_tmp_inv = inv(v_tmp,pseudo_on_fail = T)
                f2[1:n,1:n] = v_tmp_inv
                
                #tmp_m4_inv=0.25*m4(v_tmp_inv,n)
                tmp_m4_inv = 1/2*kronecker(v_tmp_inv,v_tmp_inv)
                f2[(n+1):(n+n*n),(n+1):(n+n*n)] = tmp_m4_inv
                
                # f2[1:n,1:n]=v_tmp\diag_matlab(1,size(v_tmp))
                #f2[1:n,1:n]=inv(v_tmp)
                #m4_tmp = m4(v_tmp,n)
                #f2[(n+1):(n+n*n),(n+1):(n+n*n)]=m4_tmp\diag_matlab(1,size(m4_tmp))
                #f2[(n+1):(n+n*n),(n+1):(n+n*n)]=inv(m4(v_tmp,n))
            
            ret=ret + np.matmul(np.matmul(np.transpose(f1),f2),f1)

    ret = ret/poped_db["settings"]["iFOCENumInd"]

    return {"ret": ret, "poped_db": poped_db}