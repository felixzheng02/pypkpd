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


import numpy as np
from matpy.matrix import Matrix

def diag_matlab(mat:Matrix):
    dim_mat = mat.get_shape()
    if dim_mat is not None:
        if  1 in dim_mat:
            if all(dim_mat[i] == 1 for i in range(0, len(dim_mat))) is False:
                return Matrix(np.diag(mat.get_all_data()))
    return Matrix(np.diag(mat.get_all_data()))


