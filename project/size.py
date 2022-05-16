"""
## size(obj: int/np.ndarray/pd.Dataframe/matrix) -> list
## @param obj: an object (that you want to know the dimensions). Typically a matrix.
## @return: the dimensions of the object or specific dimension you are interested in. 


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
import pandas as pd
from matpy.matrix import matrix


def size(obj): # @param dimension_index is removed
	dim_obj = None
	if type(obj) is int:
		dim_obj = None
	elif type(obj) is np.ndarray or type(obj) is pd.DataFrame or type(obj) is matrix:
		if type(obj) is matrix:
			obj = obj.get_all_data()
		dim_obj = obj.shape
		if len(dim_obj) == 1:
			dim_obj = [1, dim_obj[0]]
		else:
			dim_obj = list(dim_obj)
	if dim_obj is None:
		dim_obj = [1, 1]
	return dim_obj