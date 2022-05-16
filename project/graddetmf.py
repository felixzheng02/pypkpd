"""
## Function written to match MATLAB function gradetmf()


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.cell import cell
from project.size import size
from project.zeros import zeros
from project.mftot import mftot
from project.mf_all import mf_all
from project.get_fim_size import get_fim_size
from project.update_designinlist import update_designinlist


def graddetmf(model_switch,aX,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped_db,lndet=False,gradxt=False):
  
    n = get_fim_size(poped_db)
    m = size(ni)[0]
    if (gradxt == False):
        gdmf = np.ones(m*size(a)[1]).reshape(m, size(a)[1])
    else:
        gdmf = np.ones(m*size(xt)[1]).reshape(m, size(xt)[1])
    
    
    iParallelN = (poped_db["settings"]["parallel"]["bParallelSG"] == 1) + 1 #1 if no parallel, 2 if parallel
    
    if iParallelN == 2:
        designsin = cell(1,0)
        it = 1
    
    for p in range(0, iParallelN):
        if p == 2:
            #Execute parallel designs
            designout = designsin
            raise Exception("Parallel execution not yet implemented in PopED for R")
            #designout = execute_parallel(designsin,poped_db)
        designout = designsin


        if iParallelN == 1:
            returnArgs = mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped_db) 
            mft = returnArgs[0]
            poped_db = returnArgs[1]
        else:
            if p == 1:
                designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
            else:
                mft = designout["it"]["FIM"]
                it = it + 1
        
        
        if iParallelN == 1 or p == 2:
            if all(size(poped_db["settings"]["prior_fim"]) == size(mft)):
                mft = mft + poped_db["settings"]["prior_fim"]
            
            imft = np.linalg.inv(mft)
            if np.isinf(imft[1,1]):
                imft = zeros(size(mft))
            
        for i in range(0, m):
            if groupsize[i] == 0:
                gdmf[i, 1:ni[i]] = zeros(1, ni(i))
            else:
                if x.size != 0:
                    x_i = np.transpose(x[i,:])
                else:
                    x_i =  zeros(0,1)
                
                if a.size != 0:
                    a_i = np.transpose(a[i,:])      
                else:
                    a_i =  zeros(0,1)
                
                if iParallelN == 1:
                    returnArgs = mf_all(np.transpose(model_switch[i,1:ni[i]]), np.transpose(xt[i,1:ni[i]]),x_i,a_i,bpop,d,sigma,docc,poped_db) 
                    mf_tmp = returnArgs[0]
                    poped_db = returnArgs[1]
                else:
                    if p == 1:
                        designsin = update_designinlist(designsin, 1, ni, xt, x, a, -1, i)
                    else:
                        mf_tmp = designout["it"]["FIM"]
                        it = it + 1
                    
            mfb = groupsize[i]*mf_tmp
            
        a0 = a
        xt0 = xt
        if gradxt is False:
            nCtl = size(poped_db["design"]["a"])[1]
        else:
            nCtl = ni[i]
            for ct1 in range(0, nCtl):
                if aX[i,ct1] != 0:
                    if gradxt == False:
                        a = a0
                        a[i,ct1] = a[i,ct1] + poped_db["settings"]["hgd"]
                    else:
                        xt = xt0
                        xt[i,ct1] = xt[i,ct1] + poped_db["settings"]["hgd"]
                        
                    if iParallelN == 1:
                        returnArgs = mf_all(np.transpose(model_switch[i,1:ni[i]]), np.transpose(xt[i,1:ni[i]]), x_i, np.transpose(a[i,:]), bpop, d, sigma, docc, poped_db) 
                        mf_tmp = returnArgs[0]
                        poped_db = returnArgs[1]
                    else:
                        if p == 1:
                            designsin = update_designinlist(designsin, 1, ni, xt, x, a, -1, i)
                        else:
                            mf_tmp = designout["it"]["FIM"]
                            it = it + 1
                        
                    if iParallelN == 1 or p == 2:
                        mf_plus = groupsize[i]*mf_tmp
                        ir = (mf_plus-mfb)/poped_db["settings"]["hgd"]
                        s = 0 #Calc the tr(A^-1 * dA/dX) for some X
                        for ct2 in range(0, n):
                            s = s + np.matmul(imft[ct2,:], ir[:,ct2])
                    
                        if s == 0:  #The model doesn't depend on design variable (a or t), e.g. PD is only dependent on time and not dose, fix the a-gradient to a small value or PD with Placebo dose, fix the xt-gradient to a small value.
                            s = 1e-12
                    
                    gdmf[i, ct1] = s
                    
    if lndet == False:
        ret = gdmf*np.linalg.det(mft)
    else:
        ret = gdmf

    return {"ret": ret, "poped_db": poped_db}

