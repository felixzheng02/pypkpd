"""

Author: Caiya Zhang, Yuchen Zheng
"""

import numpy as np

def is_not_none(d: dict, k: str):
	if k in list(d.keys()):
		if d[k] is not None:
			return True
		else:
			return False
	else:
		return False

def trans(x, exp_index):
        x[exp_index] = np.exp(x[exp_index])
        return x
#trans  <- function(x) matrix(c(x[bpop_index],exp(x[d_index])),ncol=1,byrow=T

 
