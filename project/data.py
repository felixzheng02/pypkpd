import numpy as np
from matpy.matrix import Matrix
from matpy.num import Num

def data(input):
    if type(input) is Matrix:
        return input.get_data()
    elif type(input) is Num:
        return input.get_value
    elif type(input) is np.ndarray:
        if len(input.shape) == 1:
            if input.shape[0] == 0: # for shape [0,]
                reshape = [0, 1]
            else:
                reshape = [1, input.shape[0]]
            input = input.reshape(reshape)
        return input