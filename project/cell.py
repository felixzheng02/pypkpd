"""
## Create a cell array (a matrix of lists) 
 
## @param ... Dimensions for the cell array. 
## @return A list of empty lists.


## cell(mat: np.ndarray) -> int
## Returns the sum along the diagonal of the array.
## If 1-d, returns mat[0].
## If 2-d, returns np.trace(mat).
## Raises exception for other dimensions.
## @param mat: a matrix input
## @return: sum along the diagonal of the matrix, as an integer


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def cell(*argv):
	nargs = len(argv)
	dims: np.ndarray
	if nargs == 1:
		dims = np.array(argv[0])
	else:
		dims = np.array([int(i) for i in argv])
	if dims.size == 1:
		dims = np.append(dims, dims)
	if dims.size <= 1:
		raise Exception("dimensions must be of length greater than 1")
	elif any([(i < 0) for i in dims]):
		L = np.full(dims.tolist(), np.nan)
		return L
	return createCellArray(dims)


def createCellArray(dims: np.ndarray):
	L = np.full([np.prod(dims),], np.nan)
	L = np.zeros(dims.tolist())
	return L
