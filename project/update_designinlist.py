"""

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np

def update_designinlist(design,designsin,groupsize,ni,xt,x,a,it,iGroup):
  
    design["xt"] = xt
    design["a"] = a
    design["x"] = x
    design["groupsize"] = groupsize
    design["ni"] = ni
    
    if it == -1:
        design["num"] = len(designsin) + 1
    else:
        design["num"] = it
    
    design["iGroup"] = iGroup
    
    designsin = np.array([designsin, design])
    
    return designsin
