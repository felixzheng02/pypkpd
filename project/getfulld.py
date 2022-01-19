"""
## Create a full D (between subject variability) matrix given a vector of variances and covariances.
## Note, this does not test matching vector lengths.
## 
## @param variance_vector The vector of the variances.
## @param covariance_vector A vector of the covariances. Written in column major 
##   order for the lower triangular matrix.
## @return The full matrix of variances for the between subject variances
## @export
## @keywords internal
## Function written to match MATLAB function get_fulld()

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from matpy.matrix import matrix
from project.diag_matlab import diag_matlab


def getfulld(variance_vector: matrix, covariance_vector: matrix = None):
    if variance_vector.get_size() == 1:
        return variance_vector

    d = diag_matlab(variance_vector)
    if covariance_vector.get_size() > 0 and sum(covariance_vector.get_all_data() != 0) > 0:

        #d[lower.tri(d)] = covariance_vector
        d = matrix(np.transpose(d.get_all_data())) # upper.tri has wrong order, so fill lower, transpose this to upper, then fill lower again
        #d[lower.tri(d)] = covariance_vector
    
    return d 

