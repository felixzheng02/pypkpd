"""
## Function written to match MATLAB's diag function
## 
## This function tries to match the MATLAB version in handling vectors
## (matricies with one dimension equal to one), and will return a diagonal
## Matrix in these situations.
## 
## @param mat Either a vector to make into a diagonal Matrix or a Matrix you 
##   want to extract the diagonal from
## @return Either a diagonal Matrix or the diagonal of a Matrix.
## @family MATLAB
## @family Matrix_manipulation
## @example test/Test_diag_matlab.py
## @export
## @keywords internal

## Author: Caiya Zhang, Yuchen Zheng

"""


# test passed


import numpy as np
from matpy.matrix import Matrix
from project.size import size
from project.data import data
from project.length import length

def diag_matlab(mat):
    dim_mat = size(mat)
    if dim_mat is not None:
        if 1 in dim_mat:
            if np.any((np.array(dim_mat) != 1)): # mat is 1*n
                return Matrix(np.diag(data(mat).reshape(length(mat),)))
            else: # mat is 1*1
                return Matrix(np.diagonal(mat.get_data()))
    return Matrix(np.diag(data(mat)))