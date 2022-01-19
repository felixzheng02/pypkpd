"""


## Author: Caiya Zhang, Yuchen Zheng
"""


import re
import inspect
import numpy as np
from project.sfg import sfg


def find_largest_index(func_str="sfg",lab="bpop",mat=False,mat_row=True):
    if callable(func_str):
        txt = inspect.getsource(func_str)
    else: 
        txt = inspect.getsource(eval(func_str))

    txt = re.findall(('\\"[\S]*'+lab+"\\["+"[^\\,]*"), txt, re.I)
    txt = " ".join(txt)
    txt = "parameters=np.array(" + txt + ")"
    # " ".join(list(
    #     .groups()))
    ind = 0
    if len(txt) != 0 and mat is False: 
        ind = re.findall("\d", txt)
        ind = np.unique(ind)
    if len(txt) != 0 and mat and mat_row is True: 
        ind = re.findall("\d", txt)
        ind = np.unique(ind)
    if len(txt) != 0 and mat and mat_row is False: 
        ind = re.findall("\d", txt)
        ind = np.unique(ind)
    
    ind = np.array(ind).astype(float)
    if np.size(ind) == 0:
        return 0
    else:
        return np.max(ind)

# 
#  find_largest_index("sfg","bpop")
#  find_largest_index("sfg","b")
#  find_largest_index("sfg","bocc",mat=T,mat.row=T)
#  find_largest_index("sfg","x")
#  find_largest_index("sfg","a")
# 
# txt = capture.output(eval(parse(text="sfg")))
# txt = grep(paste("^[^\\#]*bpop","\\[",sep=""),txt,value=T)
# txt
# ind = gsub(paste("^[^\\#]*","bpop","\\[","(\\d+)\\].*",sep=""),"\\1",txt)
# max(as.numeric(ind))