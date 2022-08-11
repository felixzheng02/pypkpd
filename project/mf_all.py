""""
## Function written to match MATLAB function mf_all()
## Then manually adjusted to make work


## Author: Caiya Zhang, Yuchen Zheng
"""


from project.mf3 import mf3
from project.mf7 import mf7
from project.util import get_dict_value


def mf_all(model_switch, xt_ind, x, a, bpop, d, sigma, docc, pypkpd_db):

    iFIMCalculationType = get_dict_value(pypkpd_db, "settings", "iFIMCalculationType") + 1
    returnArgs = None
    if iFIMCalculationType == 1:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,pypkpd_db) #Default (with no assumption that bpop and b are uncorrelated)
    elif iFIMCalculationType == 2:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,pypkpd_db) #Reduced FIM
    elif iFIMCalculationType == 3:
        raise Exception("Not yet implemented")
    elif iFIMCalculationType == 4:
        raise Exception("Not yet implemented")
    elif iFIMCalculationType == 5:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,pypkpd_db) #Reduced FIM with derivative of SD sigma
    elif iFIMCalculationType == 6:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,pypkpd_db) #FULL FIM parameterized with A,B,C matrices & derivative of variance
    elif iFIMCalculationType == 7:
        returnArgs = mf7(model_switch,xt_ind,x,a,bpop,d,sigma,docc,pypkpd_db) #Calculate one model switch at a time, good for large matrices
    elif iFIMCalculationType == 8:
        returnArgs = mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,pypkpd_db) #Reduced FIM parameterized with A,B,C matrices & derivative of variance

    if returnArgs is None: 
        raise Exception("Unknown FIM-calculation type")
    
    return returnArgs