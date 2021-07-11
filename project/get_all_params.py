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
from project.diag_matlab import diag_matlab


def get_all_params (poped_db):
  #Return all params (in a vector all) with the specified order above
  
    bpop = poped_db["parameters"]["bpop"][:,1]
    d = poped_db["parameters"]["d"][:.1]
    docc = poped_db["parameters"]["docc"][:,1]
    covd = poped_db["parameters"]["covd"]
    sigma = diag_matlab(poped_db["parameters"]["sigma"])
    covsigma = zeros(1,(sigma.size)*(sigma.size-1)/2)
  
    k = 1

    for i in range(0, size(poped_db["parameters"]["sigma"])[0]):
        for j in range(0, size(poped_db["parameters"]["sigma"])[1]):
            if i < j:
                covsigma[k] = poped_db["parameters"]["sigma"][i,j]
                k = k + 1

  
    covdocc = poped_db["parameters"]["covdocc"]
    
    all = np.array([[bpop], [d], [np.transpose(covd)], [docc], [np.transpose(covdocc)], [sigma], [np.transpose(covsigma)]])
    return {"bpop": bpop, "d": d, "covd": covd, "docc": docc, "covdocc": covdocc, "sigma": sigma, "covsigma": covsigma, "all": all}

