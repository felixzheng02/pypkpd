"""
#' Create a full D (between subject variability) matrix given a vector of variances and covariances.
#' Note, this does not test matching vector lengths.
#' 
#' @param variance_vector The vector of the variances.
#' @param covariance_vector A vector of the covariances. Written in column major 
#'   order for the lower triangular matrix.
#' @return The full matrix of variances for the between subject variances
#' @example tests/testthat/examples_fcn_doc/examples_getfulld.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'


## Author: Andrew Hooker
"""


import numpy as np
from project.diag_matlab import diag_matlab

def getfulld(variance_vector: np.ndarray,covariance_vector: np.ndarray=None):
    if variance_vector.size == 1:
        return variance_vector

    d = diag_matlab(variance_vector)
    if covariance_vector.size > 0 and sum(covariance_vector != 0) > 0:

        #d[lower.tri(d)] = covariance_vector
        d = np.transpose(d) # upper.tri has wrong order, so fill lower, transpose this to upper, then fill lower again
        #d[lower.tri(d)] = covariance_vector
    
    return d 

