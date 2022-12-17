import numpy as np
from matpy.matrix import Matrix

def length(input) -> int:
    if type(input) is int or type(input) is np.float64:
        return 1
    elif type(input) is list:
        # only works for "completely filled list"
        result = 1
        tmp = input
        while (type(tmp[0]) is list):
            result = result * len(tmp)
            tmp = tmp[0]
        return result
    elif type(input) is np.ndarray:
        return input.size
    elif type(input) is Matrix:
        return input.get_size()