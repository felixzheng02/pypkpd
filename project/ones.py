"""
Create a Matrix of ones of size (dim1 * dim2)

@param dim1 The dimension of the Matrix (if square) or the number of rows
@param dim2 The number of columns
@return A Matrix of ones

Author: Caiya Zhang, Yuchen Zheng
"""


import path
import numpy as np
from matpy.matrix import Matrix

def ones(shape: list):
	mat = np.ones(tuple(shape))
	return Matrix(mat)