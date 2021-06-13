"""
Test if the size of a matrix (m*n) is correct

@param correct_size -> array The correct size of a matrix
@param mat The matrix being tested
@param name The name of the matrix
@return A boolean
"""


import numpy as np


def test_mat_size(correct_size: np.ndarray, mat: np.ndarray, name: str):
	if (correct_size == np.array(mat.shape)).all():
		return 1
	else:
		tmp1 = '*'.join(str(i) for i in mat.shape)
		tmp2 = '*'.join(str(i) for i in correct_size.tolist())
		raise Exception(name + "has dimensions " + tmp1 + " and should be " + tmp2)