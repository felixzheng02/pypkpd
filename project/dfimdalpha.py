"""


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from project.mftot import mftot
from project.getfulld import getfulld




def dfimdalpha(alpha, model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped_db,ha):


    bpop = bpopdescr[:,1]
    bpop[bpopdescr[:,0] != 0] = alpha[1:sum(bpopdescr[:,0] != 0)]
    d = ddescr[:,1]
    d[ddescr[:,0] != 0] = alpha[(sum(bpopdescr[:,0] != 0)+1): alpha.size]
    d = getfulld(d, covd)
    fim = mftot(model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpop, d, sigma, docc, poped_db)
    fim = fim["ret"]
    grad = np.zeros(1).reshape[size(fim)[0], size(fim)[1], alpha.size]
    for i in range(0, alpha.size):
        alpha_plus = alpha
        alpha_plus[i] = alpha_plus[i] + ha
        bpop=bpopdescr[:,1]
        bpop[bpopdescr[:,0] != 0] = alpha_plus[1:sum(bpopdescr[:,0] != 0)]
        d = ddescr[:,1]
        d[ddescr[:,0] != 0] = alpha_plus[(sum(bpopdescr[:,0] != 0) + 1): alpha_plus.size]
        d = getfulld(d, covd)
        fim_plus = mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,poped_db)
        fim_plus = fim_plus["ret"]
        if i > sum(bpopdescr[:,0]!=0):
            grad[:,:,i] = (fim_plus-fim)/ha
        else:
            #central differences for fixed effects
            alpha_minus = alpha
            alpha_minus[i] = alpha_minus[i] - ha
            bpop = bpopdescr[:,1]
            bpop[bpopdescr[:,0] != 0] = alpha_minus[1:sum(bpopdescr[:,0] != 0)]
            fim_minus = mftot(model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpop, d, sigma, docc, poped_db)
            fim_minus = fim_minus["ret"]
            grad[:,:,i] = (fim_plus-fim_minus)/(2*ha)
        
    return {"grad": grad, "fim": fim}
  