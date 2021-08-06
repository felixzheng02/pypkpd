"""
## sfg function

"""

import numpy as np
from matpy.matpy import matpy

def sfg(x,a,bpop,b,bocc):
    parameters = matpy(np.array([bpop[0]*np.exp(b[0]), 
                                bpop[1]*np.exp(b[1]),
                                bpop[2]*np.exp(b[2]),
                                bpop[3],
                                a[0]]),
                        (1,5),
                        ["CL", "V", "KA", "Favail", "DOSE"],
                        None, None)
    return parameters

