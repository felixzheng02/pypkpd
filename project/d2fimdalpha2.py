"""

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from project.mftot import mftot
from project.getfulld import getfulld


def d2fimdalpha2(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped_db,ha): 
    bpop = bpopdescr[:,1]
    bpop[bpopdescr[:,0] != 0] = alpha[1:sum(bpopdescr[:,0] != 0)]
    d = ddescr[:,1]
    d_alpha = False
    num_alpha = alpha.size
    if alpha.size > sum(bpopdescr[:,0] != 0):
        d_alpha = True
        d[ddescr[:,0] != 0] = alpha[(sum(bpopdescr[:,0] != 0) + 1):num_alpha]
    
    d = getfulld(d, covd)
    
    fim = mftot(model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpop, d, sigma, docc, poped_db)["ret"]
    
    i = 1
    hess = np.zeros([size(fim)[0],size(fim)[1],alpha.size,alpha.size])
    
    for i in range(0, alpha.size):
        alpha_plus = alpha
        alpha_plus[i] = alpha[i] + ha
        bpop = bpopdescr[:,1]
        bpop[bpopdescr[:,0] != 0] = alpha_plus[1:sum(bpopdescr[:,0] != 0)]
        if d_alpha is True:
            d = ddescr[:,1]
            d[ddescr[:,0] != 0] = alpha_plus[(sum(bpopdescr[:,0] != 0) + 1):num_alpha]
            d = getfulld(d, covd)
        
        fim_plus = mftot(model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpop, d, sigma, docc, poped_db)["ret"]
        
        for j in range(0,i):
            alpha_plus2 = alpha
            alpha_plus2[j] = alpha[j] + ha
            bpop = bpopdescr[:,1]
            bpop[bpopdescr[:,0] != 0] = alpha_plus2[1:sum(bpopdescr[:,0] != 0)]
            if d_alpha is True: 
                d = ddescr[:,1]
                d[ddescr[:,0] != 0] = alpha_plus2[(sum(bpopdescr[:,0] != 0) + 1):num_alpha]
                d = getfulld(d, covd)
            
            fim_plus2 = mftot(model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpop, d, sigma, docc, poped_db)["ret"]
            
            alpha_plus_plus = alpha
            alpha_plus_plus[i] = alpha_plus_plus[i] + ha
            alpha_plus_plus[j] = alpha_plus_plus[j] + ha
            
            bpop = bpopdescr[:,1]
            bpop[bpopdescr[:,0] != 0] = alpha_plus_plus[1:sum(bpopdescr[:,0] != 0)]
            if d_alpha is True:
                d = ddescr[:,1]
                d[ddescr[:,0] != 0] = alpha_plus_plus[(sum(bpopdescr[:,0] != 0) + 1):num_alpha]
                d = getfulld(d, covd)
            
            fim_plus_plus = mftot(model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpop, d, sigma, docc, poped_db)["ret"]
            hess[:,:,i,j] = (fim_plus_plus - fim_plus - fim_plus2 + fim)/ha^2
            hess[:,:,j,i] = hess[:,:,i,j]
        
    
    return {"hess": hess, "fim": fim} 

