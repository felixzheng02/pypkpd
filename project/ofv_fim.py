"""
#' Evaluate a criterion of the Fisher Information Matrix (FIM)
#' 
#' Compute a criterion of the FIM given the model, parameters, design and methods defined in the 
#' PopED database. 
#' 
#' @param fmf The FIM
#' @param poped_db A poped database
#' @param ofv_calc_type  OFV calculation type for FIM
#' \itemize{ 
#' \item 1 = "D-optimality". Determinant of the FIM: det(FIM)
#' \item 2 = "A-optimality".  Inverse of the sum of the expected parameter variances: 
#' 1/trace_matrix(inv(FIM)) 
#' \item 4 = "lnD-optimality".  Natural logarithm of the determinant of the FIM: log(det(FIM)) 
#' \item 6 = "Ds-optimality". Ratio of the Determinant of the FIM and the Determinant of the uninteresting
#' rows and columns of the FIM: det(FIM)/det(FIM_u)
#' \item 7 = Inverse of the sum of the expected parameter RSE: 1/sum(get_rse(FIM,poped_db,use_percent=FALSE))
#' }
#' @param use_log Should the criterion be in the log domain?
#' @inheritParams RS_opt
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' 
#' @return The specified criterion value.
#' 
#' @family FIM
#' @family evaluate_FIM
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_ofv_fim.R
#' @export
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from project.get_cv import get_rse
from project.trace_matrix import trace_matrix


def ofv_fim(fmf,poped_db,*args):

    ofv_calc_type = poped_db["settings"]["ofv_calc_type"]
    ds_index=poped_db["parameters"]["ds_index"]
    use_log = True

    #Input: the FIM
    #Return the single value that should be maximized
    
    ## create ds_index vector if not already done
    if type(ds_index) is np.ndarray:
        ds_index = np.array([ds_index,1,ds_index.size])
    
    ofv_value = 0
    
    if len(poped_db["settings"]["prior_fim"]) != 0 and all(size(poped_db["settings"]["prior_fim"])==size(fmf)):
        fmf = fmf + poped_db["settings"]["prior_fim"]
    
    
    if ofv_calc_type==1:    #determinant of FIM
        if len(fmf) == 0:
            ofv_value = 0
        else:
            ofv_value = np.linalg.det(fmf)
        return( ofv_value ) 
    
    if ofv_calc_type == 2:  #trace of the inverse of FIM
        imf = np.linalg.inv(fmf)
        ofv_value = trace_matrix(imf)
        ofv_value = 1/ofv_value #Make it a max-problem
        return( ofv_value ) 
    
    
    if ofv_calc_type == 3:  #S-Optimal Design
        raise Exception("S-optimal design not implemented yet!")
    
    
    if ofv_calc_type==4:   #log determinant of FIM
        #ofv_value = sum(log(svd(fmf)))
        det_fim = np.linalg.det(fmf)
        if det_fim < 0: 
            det_fim = 0
        ofv_value = np.log(det_fim)
    
    
    if ofv_calc_type == 5:  # C-optimal design
        raise Exception("C-optimal design is not implemented yet!")
    

    if ofv_calc_type == 6:  #Ds-optimal design
        tmp = fmf
        ##col(ds_index) returns col index of each elem
        ##col(ds_index)[ds_index==1] 
        tmp = tmp[np.where(ds_index==1),:]
        tmp = tmp[:,np.where(ds_index==1)]
        if use_log is True:
            ofv_value = np.log(np.linalg.det(fmf))-np.log(np.linalg.det(tmp))
        else:
            ofv_value = np.linalg.det(fmf)/np.linalg.det(tmp)  
        
    if ofv_calc_type == 7:  #sum of CV
        if sum(sum(fmf)) != 0 and np.isnan(sum(sum(fmf))) is False:
            #imf = inv(fmf)
            #returnArgs <-  get_cv(diag_matlab(imf),poped_db["parameters"]bpop,poped_db["parameters"]d,poped_db["parameters"]docc,poped_db["parameters"]sigma,poped_db) 
            #params <- returnArgs[[1]]
            #params_cvs <- returnArgs[[2]]
            params_cvs = get_rse(fmf,poped_db,use_percent=False)
            if np.isnan(sum(params_cvs)) is True:
                ofv_value = 0
            else:
                ofv_value=1/sum(np.abs(params_cvs))
        else:
            ofv_value = 0
        return ofv_value 
    
    return ofv_value