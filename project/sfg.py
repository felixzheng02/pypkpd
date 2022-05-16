"""
# sfg function

"""

import numpy as np
from matpy.num import num
from matpy.matrix import matrix

def sfg(x,a,bpop,b,bocc):
    parameters = matrix(np.array([bpop.get_by_index(0)*np.exp(b.get_by_index(0)), 
                                bpop.get_by_index(1)*np.exp(b.get_by_index(1)),
                                bpop.get_by_index(2)*np.exp(b.get_by_index(2)),
                                bpop.get_by_index(3),
                                a.get_by_index(0)]),
                        (1,5),
                        ["CL", "V", "KA", "Favail", "DOSE"],
                        None, None)
    return parameters