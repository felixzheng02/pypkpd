import path
import numpy as np
from matpy.matrix import Matrix


def set_data_by_flattened_index(mat: Matrix or np.ndarray, flattened_index: int, new_data):
    """
    makes use of np.unravel_index
    to set one data based on the flattened index.
    """
    if type(mat) is Matrix:
        mat.set_one_data(new_data, index=list(np.unravel_index(flattened_index, mat.get_shape())))
    elif type(mat) is np.ndarray:
        mat[np.unravel_index(flattened_index, mat.shape)] = new_data
    return mat