"""
## Function written to match MATLAB function write_matric()


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np

def write_Matrix(f, x: np.ndarray):
	for i in range(0, x.shape[0]):
		print("%6e" % x[i, :], file = f)
		print("\n", file = f)
	return
