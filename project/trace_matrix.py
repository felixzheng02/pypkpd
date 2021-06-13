"""
Function written to match MATLAB function

Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def trace_matrix (mat):

    if len(mat.shape) == 1:
        return mat[0]
    elif len(mat.shape) == 2:
        return(np.trace(mat))
    else:
        raise Exception("Error: mat needs to be an array of 1 or 2 dimension")