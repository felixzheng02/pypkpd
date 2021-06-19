"""
@param obj An object you want to know the various dimensions of. Typically a matrix.
@param dimension_index Which dimension you are interested in.
@return The dimensions of the object or specific dimension you are interested in. 
Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
import pandas as pd


def size(obj, dimension_index=None):
	dim_obj = None
	if type(obj) is int:
		dim_obj = None
	elif type(obj) is np.ndarray or type(obj) is pd.DataFrame:
		dim_obj = obj.shape
		if len(dim_obj) == 1:
			dim_obj = [1, dim_obj[0]]
		else:
			dim_obj = list(dim_obj)
	if dim_obj is None:
		dim_obj = [1, 1]
	if dimension_index is None:
		return dim_obj
	else:
		return dim_obj[dimension_index - 1]