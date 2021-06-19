import numpy as np
from project.zeros import zeros
from project.size import size


l1 = size(obj=np.array([2,3,4,5,6]))

l2 = size(10)

l3 = size(zeros(4,7))

print(l1)
print(l2)
print(l3)