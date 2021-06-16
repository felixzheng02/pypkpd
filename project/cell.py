"""
Create a cell array (a matrix of lists) 
 
@param ... Dimensions for the cell array. 
@return A list of empty lists.

Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def cell(*argv):
	nargs = len(argv)
	dims = np.array([])
	if nargs == 1:
		dims = np.array([argv[0]])
	else:
		dims = np.array([int(i) for i in argv])
	if len(dims) == 1:
		dims = np.array(dims.tolist().append(dims[0]))
	if len(dims) <= 1:
		raise Exception("dimensions must be of length greater than 1")
	elif any(i < 0 for i in dims):
		L = np.full(dims.tolist(), np.nan)
		return L
	return createCellArray(dims)


def createCellArray(dims: np.ndarray):
	L = np.full([np.prod(dims),], np.nan)
	L = np.zeros(dims.tolist())
	return L