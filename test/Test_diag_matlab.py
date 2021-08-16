
import numpy as np
from matpy.matrix import matrix
from project.diag_matlab import diag_matlab


a = matrix(np.array([1,2,3]))
b = matrix(np.array([1,2,3,4,5,6,7,8,9]),(3,3))
print(diag_matlab(a).get_all_data())
print(diag_matlab(b).get_all_data())