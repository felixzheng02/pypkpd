"""

Author: Caiya Zhang, Yuchen Zheng
"""



def find_largest_index(func_str="sfg",lab="bpop",mat=False,mat_row=True):
    if callable(func_str):
        txt = capture_output(func_str)
    else:
        txt = capture_output(eval(parse(text=func_str)))
    
    txt = grep(paste("^[^\\#]*",lab,"\\[",sep=""),txt,value=T)
    ind = 0
    if len(txt) != 0 and mat is False: 
        ind = gsub(paste("^[^\\#]*",lab,"\\[\\s*(\\d+)\\s*\\].*",sep=""),"\\1",txt)
    if len(txt) != 0 and mat and mat_row is True: 
        ind = gsub(paste("^[^\\#]*",lab,"\\[\\s*(\\d+)\\s*,.*?\\].*$",sep=""),"\\1",txt)
    if len(txt) != 0 and mat and mat_row is False: 
        ind = gsub(paste("^[^\\#]*",lab,"\\[.*?,\\s*(\\d+)\\s*\\].*",sep=""),"\\1",txt)
    
    max(float(ind))

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
