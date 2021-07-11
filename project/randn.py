"""
#' Function written to match MATLAB's randn function
#' 
#' Generate random samples from a standardized normal distribution and return in matrix form.
#' 
#' @param dim1 The dimension of the matrix (if square), otherwise the number of rows.
#' @param dim2 The number of columns, if different from the number of rows.
#' 
#' @return Matrix of random generated samples.
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_randn.R
#' @export
#' @keywords internal

## Function written to match MATLAB function

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np

def randn (dim1, dim2=None):
    if dim2 is None:
        dim2 = dim1
    tmp = np.random.normal(dim1*dim2)
    mat = np.ndarray(tmp).reshape([dim1,dim2])
    #mat = matrix(tmp,dim1,dim2)
    return mat