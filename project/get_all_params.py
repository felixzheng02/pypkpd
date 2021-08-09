"""
#' Extract all model parameters from the PopED database.
#' 
#' @param poped_db A PopED database.
#' @return A list containing:
#' \item{bpop}{A vector of fixed effect parameter values.}
#' \item{d}{A vector of between subject variability parameters}
#' \item{covd}{A vector of the covariances of the between subject variability parameters.  Row major format of the lower triangular portion of the D (OMEGA) matrix}
#' \item{docc}{A vector of the between occasion variability (BOV) terms in the model}
#' \item{covdocc}{A vector of the covariances between the BOV terms.  Row major of the lower triangular portion of the BOV matrix. }
#' \item{sigma}{A vector of the residual unexplained variances (RUV)}
#' \item{covsigma}{A vector of the covariances between the RUV terms}
#' \item{all}{A vector with all of the above, in the order of this list.}
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_get_all_params.R
#' @export
#' @keywords internal

## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from project.zeros import zeros
from matpy.matrix import matrix
from project.diag_matlab import diag_matlab


def get_all_params (poped_db):
    #Return all params (in a vector all) with the specified order above

    #type: matrix to ndarray
    bpop = poped_db["parameters"]["bpop"].get_data()[:,1]
    d = poped_db["parameters"]["d"].get_data()[:,1]
    docc = poped_db["parameters"]["docc"].get_data()[:,1]
    covd = poped_db["parameters"]["covd"].get_data()
    covdocc = poped_db["parameters"]["covdocc"].get_data()
    sigma = diag_matlab(poped_db["parameters"]["sigma"]).get_data()
    covsigma = zeros(1,(sigma.size)*(sigma.size-1)/2).get_data()

    k = 1

    for i in range(0, size(poped_db["parameters"]["sigma"])[0]):
        for j in range(0, size(poped_db["parameters"]["sigma"])[1]):
            if i < j:
                covsigma[k] = poped_db["parameters"]["sigma"][i,j]
                k = k + 1
    
    all = matrix(np.array([[bpop], [d], [np.transpose(covd)], [docc], [np.transpose(covdocc)], [sigma], [np.transpose(covsigma)]]))
    
    all_mat = matrix(np.array([[bpop], [d], [covd], [docc], [covdocc], [sigma], [covsigma], [all]]))
    all_mat.set_datanam(["bpop", "d", "covd", "docc", "covdocc", "sigma", "covsigma", "all"])
    


    return all_mat