"""
Create a cell array (a matrix of lists) 
 
@param ... Dimensions for the cell array. 
@return A list of empty lists.
"""


import project.all_modules as am


def cell(*argv):
	nargs = len(argv)
	dims = am.np.array([])
	if nargs == 1:
		dims = am.np.array([argv[0]])
	else:
		dims = am.np.array([int(i) for i in argv])
	if len(dims) == 1:
		dims = am.np.array(dims.tolist().append(dims[0]))
	if len(dims) <= 1:
		raise Exception("dimensions must be of length greater than 1")
	elif any(i < 0 for i in dims):
		L = am.np.full(dims.tolist(), am.np.nan)
		return L
	return createCellArray(dims)


def createCellArray(dims: am.np.ndarray):
	L = am.np.full([am.np.prod(dims),], am.np.nan)
	L = am.np.zeros(dims.tolist())
	return L
