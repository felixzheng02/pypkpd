"""
Author: Caiya Zhang, Yuchen Zheng
"""


import re
import inspect


def find_largest_index(func_str="sfg",lab="bpop",mat=False,mat_row=True):
    if callable(func_str):
        txt = inspect.getsource(func_str)
    else:
        txt = inspect.getsource(eval(func_str + "()"))

    txt = re.search("^[^\\#]*" + lab + "\\[", txt)
    txt = " ".join(list(txt.groups()))
    # " ".join(list(
    #     .groups()))
    ind = 0
    if len(txt) != 0 and mat is False: 
        ind = re.sub("^[^\\#]*" + lab + "\\[\s*(\d+)\s*\\].*", "\\1", txt)
    if len(txt) != 0 and mat and mat_row is True: 
        ind = re.sub("^[^\\#]*" + lab + "\\[\s*(\d+)\s*,.*?\\].*$", "\\1", txt)
    if len(txt) != 0 and mat and mat_row is False: 
        ind = re.sub("^[^\\#]*" + lab + "\\[.*?,\s*(\d+)\\s*\\].*", "\\1", txt)
    
    return (float(ind))

# 
#  find.largest.index("sfg","bpop")
#  find.largest.index("sfg","b")
#  find.largest.index("sfg","bocc",mat=T,mat.row=T)
#  find.largest.index("sfg","x")
#  find.largest.index("sfg","a")
# 
# txt = capture.output(eval(parse(text="sfg")))
# txt = grep(paste("^[^\\#]*bpop","\\[",sep=""),txt,value=T)
# txt
# ind = gsub(paste("^[^\\#]*","bpop","\\[","(\\d+)\\].*",sep=""),"\\1",txt)
# max(as.numeric(ind))