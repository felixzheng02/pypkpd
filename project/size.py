"""
@param obj An object you want to know the various dimensions of. Typically a matrix.
@param dimension_index Which dimension you are interested in.
@return The dimensions of the object or specific dimension you are interested in. 
Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
import pandas as pd


def size(obj): # @param dimension_index is removed
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
	return dim_obj