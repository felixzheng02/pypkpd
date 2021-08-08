import numpy as np
from project.size import size
from project.zeros import zeros
from matpy.matrix import matrix


l1 = size(obj=matrix(np.array([2,3,4,5,6])))

l2 = size(10)

l3 = size(zeros(4,7))

print(l1)
print(l2)
print(l3)