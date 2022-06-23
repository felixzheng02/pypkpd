"""
zeros.R
Create a Matrix of zeros of size (dim1 x dim2)
@param dim1 The dimension of the Matrix (if square) or the number of rows
@param dim2 The number of columns
@return A Matrix of zeros
Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from matpy.matrix import Matrix


def zeros(dim1, dim2 = None):
	if type(dim1) is int and type(dim2) is int:
		return Matrix(np.zeros([dim1, dim2]))
	else:
		return Matrix(np.zeros([int(dim1), int(dim2)]))
	