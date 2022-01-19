"""

## can be deleted!!
## Function written to match MATLAB function numel()


## Author: Caiya Zhang, Yuchen Zheng
"""

import numpy as np
from project.size import size

def numel(arr):
    return np.prod(size(arr))
