"""
# sfg function

"""

import path
import numpy as np
from matpy.num import Num
from matpy.matrix import Matrix

class Sfg:
    def __init__(self) -> None:
        pass

    def sfg_1(x, a, bpop, b, bocc):
        if type(bpop) == int:
            bpop = [bpop] * 4
        if type(b) == int:
            b = [b] * 3
        parameters = Matrix(np.array([bpop[0]*np.exp(b[0]), 
                                    bpop[1]*np.exp(b[1]),
                                    bpop[2]*np.exp(b[2]),
                                    bpop[3],
                                    a]),
                            [1,5],
                            None,
                            [None, ["CL", "V", "KA", "Favail", "DOSE"]])
        return parameters
    
    def sfg_2(x, a, bpop, b, bocc):
        parameters = Matrix(np.array([bpop[0]*np.exp(b[0]), 
                                    bpop[1]*np.exp(b[1]),
                                    bpop[2]*np.exp(b[2]),
                                    bpop[3],
                                    a[0]]),
                            [1,5],
                            None,
                            [None, ["CL", "V", "KA", "Favail", "DOSE"]])
        return parameters