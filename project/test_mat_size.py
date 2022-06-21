"""
## test_mat_size(correct_size: list, input_size: list, name: str) -> int/exception
## Test if the size of a matrix (m*n) is correct.
## Return 1 if correct.
## Raise exception if incorrect.
## @param correct_size: The correct size of a matrix
## @param input_size: The input size being tested
## @param name: The name of the matrix being tested
## @return int


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def test_mat_size(correct_size: list, input_size: list, name: str):
	if correct_size == input_size:
		return 1
	else:
		tmp1 = '*'.join(str(i) for i in input_size)
		tmp2 = '*'.join(str(i) for i in correct_size)
		raise Exception(name + " has dimensions " + tmp1 + " and should be " + tmp2)