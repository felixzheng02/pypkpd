import path
import numpy as np
from matpy.matrix import Matrix
from project.getfulld import getfulld


a = getfulld(Matrix([1, 2, 3]))
print(a.get_data())
b = getfulld(Matrix([1, 2, 3]), Matrix([7, 6, 5]))
print(b.get_data())