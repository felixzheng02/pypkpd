"""
## test_mat_size(correct_size: np.ndarray, mat: np.ndarray, name: str) -> int/exception
## Test if the size of a matrix (m*n) is correct.
## Return 1 if correct.
## Raise exception if incorrect.
## @param correct_size -> array The correct size of a matrix
## @param mat The matrix being tested
## @param name The name of the matrix
## @return A boolean


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def test_mat_size(correct_size: np.ndarray, mat: np.ndarray, name: str):
	shape = list(mat.shape)
	if len(shape) == 1:
		shape = [1, shape[0]]
	if type(correct_size == np.array(shape)) is bool:
		if correct_size == np.array(shape):
			return 1
	elif (correct_size == np.array(shape)).all():
		return 1
	else:
		tmp1 = '*'.join(str(i) for i in shape)
		tmp2 = '*'.join(str(i) for i in correct_size.tolist())
		raise Exception(name + " has dimensions " + tmp1 + " and should be " + tmp2)