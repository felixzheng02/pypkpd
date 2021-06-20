"""
#' Model linearization with respect to epsilon and eta.
#' 
#' The function performs a linearization of the model with respect to the residual variability and 
#' then the between subject variability.
#' Derivative of model w.r.t. eps then eta, evaluated at eps=0 and b=b_ind.
#' 
#' @inheritParams mftot
#' @inheritParams LinMatrixH
#' @param NumEPS The number of eps() terms in the model.
#' 
#' @return A matrix of size (samples per individual x (number of sigma x number of omega)) 
#'  
#' @family FIM
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_LinMatrixLH.R
#' @export
#' @keywords internal
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work


## Author: Caiya Zhang, Yuchen Zheng
"""


from project.size import size
from project.zeros import zeros
from project.feval import feval
from project.LinMatrixH import LinMatrixH


def LinMatrixLH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,NumEPS,poped_db):
    #----------Model linearization with respect to epsilon.
    #
    # size of return is (samples per individual x (number of sigma x number of omega)) 
    #
    # derivative of model w$r.t. sigma then eta, eval at e=0 and eta
    #
    y = zeros(size(xt_ind,1), poped_db["parameters"]["NumRanEff"]*NumEPS)
    if poped_db["settings"]["iApproximationMethod"] == 0 or poped_db["settings"]["iApproximationMethod"] == 1: #No interaction
        return [y, poped_db]
    if poped_db["parameters"]["NumRanEff"] == 0:
        return [y, poped_db]
    if poped_db["settings"]["hle_switch"] == 20:
        raise Exception("Analytic derivative with interaction is not yet available!")
    
    for i in range(0, poped_db["parameters"]["NumRanEff"]):
        b_ind_plus = b_ind
        b_ind_minus=b_ind
        b_ind_plus[i] = b_ind_plus[i] + poped_db["settings"]["hle"]
        b_ind_minus[i]= b_ind_minus[i] - poped_db["settings"]["hle"]
        returnArgs =  LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,poped_db) 
        lin_plus = returnArgs[[0]]
        poped_db = returnArgs[[1]]
        returnArgs =  LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,poped_db) 
        lin_minus = returnArgs[[0]]
        poped_db = returnArgs[[1]]
        temp = (lin_plus-lin_minus)/(2*poped_db["settings"]["hle"])
###!!!!!!!!!!!
        y[,((i-1)*NumEPS+1):(i*NumEPS)]=temp[,1:NumEPS,drop=F]

    return [y, poped_db]
    

    #Helper function to get the hessian for the AD derivative
def new_ferror_file(model_switch,deriv_vec,xt_ind,x,a,bpop,bocc_ind,poped_db):
###!!!!!!!!!!!
    fg0=feval(poped_db["model"]["fg_pointer"],x,a,bpop,deriv_vec(0:poped_db["parameters"]["NumRanEff"]),bocc_ind) #Interaction
    returnArgs = feval(poped_db["model"]["ferror_pointer"],model_switch,xt_ind,fg0,deriv_vec(poped_db["parameters"]["NumRanEff"]+1:length(deriv_vec)),poped_db) 
    f_error = returnArgs[[0]]
    poped_db = returnArgs[[1]]
    return [f_error, poped_db]
