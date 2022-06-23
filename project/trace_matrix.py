"""
trace_Matrix(mat: np.ndarray) -> int
Returns the sum along the diagonal of the array.
If 1-d, returns mat[0].
If 2-d, returns np.trace(mat).
Raises exception for other dimensions.
@param mat: a Matrix input
@return: sum along the diagonal of the Matrix, as an integer

Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def trace_Matrix (mat):

    if len(mat.shape) == 1:
        return mat[0]
    elif len(mat.shape) == 2:
        return(np.trace(mat))
    else:
        raise Exception("Error: mat needs to be an array of 1 or 2 dimension")