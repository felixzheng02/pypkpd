"""
#' Function written to match MATLAB's diag function
#' 
#' This function tries to match the MATLAB version in handling vectors
#' (matricies with one dimension equal to one), and will return a diagonal
#' matrix in these situations.
#' 
#' @param mat Either a vector to make into a diagonal matrix or a matrix you 
#'   want to extract the diagonal from
#' @return Either a diagonal matrix or the diagonal of a matrix.
#' @family MATLAB
#' @family matrix_manipulation
#' @example test/Test_diag_matlab.py
#' @export
#' @keywords internal

Author: Caiya Zhang, Yuchen Zheng

"""


import numpy as np
from matpy.matrix import matrix

def diag_matlab(mat:matrix):
    dim_mat = mat.get_shape()
    if dim_mat is not None:
        if  1 in dim_mat:
            if all(dim_mat[i] == 1 for i in range(0, len(dim_mat))) is False:
                return matrix(np.diag(mat.get_data()))
    return matrix(np.diag(mat.get_data()))


