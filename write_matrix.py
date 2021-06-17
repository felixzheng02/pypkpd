"""
Function written to match MATLAB function

Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.fprintf import fprintf

def write_matrix(f, x: np.ndarray):
	for i in range(0, x.shape[0]):
		fprintf(f,"%6e",x[i, :])
        fprintf(f,"\n")
	return
