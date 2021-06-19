"""
ones.R

Create a matrix of ones of size (dim1 * dim2)

@param dim1 The dimension of the matrix (if square) or the number of rows
@param dim2 The number of columns
@return A matrix of ones

Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def ones(dim1, dim2 = None):
	if dim2 is None:
		if dim1.size == 2:
			tmp = dim1
			dim1 = tmp[0]
			dim2 = tmp[1]
		elif dim1.size == 1:
			dim2 = dim1
		else: # stop("first argument can only have one or two values")
			print("first argument can only have one or two values")
			return
	mat = np.ones([dim1, dim2])
	return mat