"""


## Function written to match MATLAB function
## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size

def numel(arr):
    return np.prod(size(arr))
