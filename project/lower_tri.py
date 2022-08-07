"""
Returns a matrix of logicals the same size of a given matrix 
with entries True in the lower triangle.
"""


from matpy.matrix import Matrix
import numpy as np


def lower_tri(mat: Matrix):
    return np.tril(np.ones(mat.get_shape()), k=-1) == 1