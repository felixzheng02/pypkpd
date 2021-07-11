""""
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work


## Author: Caiya Zhang, Yuchen Zheng
"""


from project.mf3 import mf3
from project.mf7 import mf7



def mf_all(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped_db):

    if poped_db["settings"]["iFIMCalculationType"] + 1 == 1:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped_db) #Default (with no assumption that bpop and b are uncorrelated)
    elif poped_db["settings"]["iFIMCalculationType"] + 1 == 2:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped_db) #Reduced FIM
    elif poped_db["settings"]["iFIMCalculationType"] + 1 == 3:
        raise Exception("Not yet implemented")
    elif poped_db["settings"]["iFIMCalculationType"] + 1 == 4:
        raise Exception("Not yet implemented")
    elif poped_db["settings"]["iFIMCalculationType"] + 1 == 5:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped_db) #Reduced FIM with derivative of SD sigma
    elif poped_db["settings"]["iFIMCalculationType"] + 1 == 6:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped_db) #FULL FIM parameterized with A,B,C matrices & derivative of variance
    elif poped_db["settings"]["iFIMCalculationType"] + 1 == 7:
        returnArgs = mf7(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped_db) #Calculate one model switch at a time, good for large matrices
    elif poped_db["settings"]["iFIMCalculationType"] + 1 == 8:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped_db) #Reduced FIM parameterized with A,B,C matrices & derivative of variance

    
    if returnArgs is None: 
        raise Exception("Unknown FIM-calculation type")
    ret = returnArgs[[0]]
    poped_db = returnArgs[[1]]
    
    return {"ret": ret, "poped_db": poped_db} 



