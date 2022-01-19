"""
## Function written to match MATLAB function graddetmf_ext()


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.cell import cell
from project.size import size
from project.zeros import zeros
from project.mftot import mftot
from matpy.matrix import matrix
from project.get_fim_size import get_fim_size
from project.update_designinlist import update_designinlist


def graddetmf_ext(model_switch,aX,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped_db,lndet=False,gradxt=False):
  
    n = get_fim_size(poped_db)
    m = size(ni)[0]
    if gradxt is False:
        gdmf = matrix(1, (m, size(a)[1]))
        G_X = poped_db["design_space"]["G_a"]
    else:
        gdmf = matrix(1, (m, size(xt)[1]))
        G_X = poped_db["design_space"]["G_xt"]
    
    
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
        if iParallelN == 1:
            returnArgs = mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped_db) 
            mft = returnArgs[0]
            poped_db = returnArgs[1]
        else:
            if p == 1:
                designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
            else:
                designout = designsin
                mft = designout[it]["FIM"]
                it = it + 1


        if iParallelN == 1 or p == 2:
            #If we have a prior
            if all(size(poped_db["settings"]["prior_fim"]) == size(mft)):
                mft = mft + poped_db["settings"]["prior_fim"]
            
            imft = np.linalg.inv(mft)
            if imft[1,1] is np.inf:
                imft = zeros(size(mft))
            
        
        a0 = a
        xt0 = xt
        for k in range(0, max(max(max(G_X)), 0)):
            inters = (G_X == k)
            if sum(sum(inters)) !=0 : #If we have a covariate or time-point defined here (accord. to G_X)
                if gradxt is False:
                    a = a0 + poped_db["settings"]["hgd"]*inters
                else:
                    xt = xt0 + poped_db["settings"]["hgd"]*inters
                    
                if iParallelN == 1:
                    returnArgs =  mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped_db) 
                    mf_plus = returnArgs[0]
                    poped_db = returnArgs[1]
                else:
                    if p == 1:
                        designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
                    else:
                        mf_plus = designout[it]["FIM"]
                        it = it + 1
                    
                    if iParallelN == 1 or p == 2:
                        #If we have a prior
                        if all(size(poped_db["settings"]["prior_fim"]) == size(mft)):
                            mf_plus = mf_plus + poped_db["settings"]["prior_fim"]
                        
                        ir = (mf_plus-mft)/poped_db["settings"]["hgd"]
                        
                        s = 0 #Calc the tr(A^-1 * dA/dX) for some X
                        for ct2 in range(0, n):
                            s = s + np.matmul(imft[ct2,:], ir[:,ct2])
                        
                        
                        if s == 0: #The model doesn't depend on a or xt, e$g. PD is only dependent on time and not dose, fix the a-gradient to a small value or PD with Placebo dose, fix the xt-gradient to a small value
                            s = 1e-12
                        
                        gdmf[inters==1 and aX!=0] = s
                        # for(i in 1:size(a,1)){
                        #   for(j in 1:size(a,2)){
                        #     if((inters[i,j]==1 && aX[i,j]!=0)){
                        #       gdmf[i,j]=s
                        #     }
                        #   }
                        # }
                        
    
    if lndet is False:
        ret = gdmf*np.det(mft)
    else:
        ret = gdmf
    
    return {"ret": ret, "poped_db": poped_db}



