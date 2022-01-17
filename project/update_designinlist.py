"""
update_designinlist(design: dict,
                    designsin,
               |----groupsize,
               |    ni,
all            |    xt,
design         |    x,
keys           |    a,
               |    it,
               |----iGroup
                    ) -> np.ndarray
@param design
@param designsin
@params groupsize, ni, xt, x, a, it, iGroup: these are all key values in design dict
@return: "designsin" as an np.ndarray

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