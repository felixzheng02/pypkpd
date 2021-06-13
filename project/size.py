"""
@param obj An object you want to know the various dimensions of. Typically a matrix.
@param dimension_index Which dimension you are interested in.
@return The dimensions of the object or specific dimension you are interested in. 

Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def size(obj, dimension_index=None):
	dim_obj = None
	if type(obj) == int:
		dim_obj = None
	elif type(obj) == np.ndarray:
		dim_obj = obj.shape
	if dim_obj == None:
		dim_obj = np.array([1, 1])
	if dimension_index == None:
		return dim_obj
	return dim_obj[dimension_index + 1]