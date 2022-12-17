"""
Author: Caiya Zhang, Yuchen Zheng
"""


import path
import numpy as np
from matpy.matrix import Matrix


def size(input) -> list:
    if type(input) is list:
        # only works for "completely filled list"
        result = [1]
        tmp = input
        while (type(tmp[0]) is list):
            result.append(len(tmp))
            tmp = tmp[0]
        return result
    elif type(input) is np.ndarray:
        return input.shape
    elif type(input) is Matrix:
        return input.get_shape()
    elif input is None:
        return [1, 0]
    else: # if it is int or float
        return [1, 1]