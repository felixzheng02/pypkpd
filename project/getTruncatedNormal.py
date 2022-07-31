"""
##  Generate a random sample from a truncated normal distribution.
##  
## @param mean the mean of the normal distribution
## @param variance The variance of the normal distribution
##    
## @return A random sample from the specified truncated normal distribution
## 
## @export
## @keywords internal


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np

def getTruncatedNormal(mean, variance):
    while True:
        n = mean + np.random.randn(1, 1) * np.sqrt(variance)
        if np.sign(n) == np.sign(mean): # n is 1*1, only 1 element
            break
    return n