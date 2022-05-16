"""
## subspace_min(x, l, u, x_cp, d, G) -> np.ndarray (?)
## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np


def subspace_min (x, l, u, x_cp, d, G):
    n =  len(x)
    Z = np.diag(n)
    fixed = (x_cp <= l+1e-8) or (x_cp >= u-1e8)
    if all(fixed) is True:
        x = x_cp
        return x
    
    Z = Z[: , fixed is False]
    rgc = np.matmul(np.transpose(Z), (d + np.matmul(G, (x_cp-x))))
    rB = np.matmul(np.matmul(np.transpose(Z), G), Z)
    d[fixed is False] = np.linalg.solve(rB, rgc)
    d[fixed is False] = -d[fixed is False]
    alpha = 1
    temp1 = alpha
    for i in np.where(fixed is False):
        dk = d[i]
        if dk < 0:
            temp2 = l[i] - x_cp[i]
            if temp2 >= 0:
                temp1 = 0
            else:
                if dk*alpha < temp2: 
                    temp1 = temp2/dk
            
        else:
            temp2 = u[i] - x_cp[i]
            if temp1 <= 0:
                temp1 = 0
            else:
                if dk*alpha > temp2: 
                    temp1 = temp2/dk
            
        alpha = min(temp1, alpha)
    
    x = x_cp + Z*alpha*d[fixed is False]
    return x