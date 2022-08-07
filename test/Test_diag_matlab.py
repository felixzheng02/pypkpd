import path
import numpy as np
from matpy.matrix import Matrix
from project.diag_matlab import diag_matlab


a = diag_matlab(Matrix(3))
print(a.get_data())
b = diag_matlab(Matrix([1, 2, 3]))
print(b.get_data())
c = diag_matlab(Matrix([1, 2, 3], shape=[3, 1]))
print(c.get_data())
d = diag_matlab(Matrix([1, 2, 3, 4], shape=[2, 2]))
print(d.get_data())