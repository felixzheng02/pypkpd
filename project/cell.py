"""
@param Dimensions
@return np.ndarray filled with np.nan of that dimension

Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def cell(shape: list, fill_value = np.nan):
	return np.full(shape=shape, fill_value=fill_value)