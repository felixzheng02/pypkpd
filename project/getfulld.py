"""
## Create a full D (between subject variability) Matrix given a vector of variances and covariances.
## Note, this does not test matching vector lengths.
## 
## @param variance_vector The vector of the variances.
## @param covariance_vector A vector of the covariances. Written in column major 
##   order for the lower triangular Matrix.
## @return The full Matrix of variances for the between subject variances
## @export
## @keywords internal
## Function written to match MATLAB function get_fulld()

## Author: Caiya Zhang, Yuchen Zheng
"""


# test passed


import numpy as np
from matpy.matrix import Matrix
from project.diag_matlab import diag_matlab
from project.lower_tri import lower_tri


def getfulld(variance_vector: Matrix, covariance_vector: Matrix = None):
    if variance_vector.get_size() == 1:
        return variance_vector

    d = diag_matlab(variance_vector)
    if covariance_vector is not None:
        if covariance_vector.get_size() > 0 and np.sum(covariance_vector.get_data() != 0) > 0:

            d.data[lower_tri(d)] = covariance_vector.get_data().reshape(covariance_vector.get_size(),)
            d = Matrix(np.transpose(d.get_data())) # upper.tri has wrong order, so fill lower, transpose this to upper, then fill lower again
            d.data[lower_tri(d)] = covariance_vector.get_data().reshape(covariance_vector.get_size(),)
    
    return d 