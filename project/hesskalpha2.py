"""

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from project.zeros import zeros
from project.dfimdalpha import dfimdalpha
from project.d2fimdalpha2 import d2fimdalpha2
from project.trace_Matrix import trace_Matrix
from project.log_prior_pdf import log_prior_pdf



def hesskalpha2(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped_db,ha,Engine):
    #D2KALPHA2 calculates the hessian of k with respect to alpha
    #   Detailed explanation goes here
    returnArgs = log_prior_pdf(alpha, bpopdescr, ddescr, return_gradient=True, return_hessian=True) 
    p = returnArgs[0]
    gradp = returnArgs[1]
    hessp = returnArgs[2]
    #get dF/dAlpha and fim
    returnArgs = dfimdalpha(alpha, model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpopdescr, ddescr, covd, sigma, docc, poped_db, ha) 
    d_fim = returnArgs[0]
    fim = returnArgs[1]
    ifim = np.linalg.inv(fim)
    tigi = zeros(size(d_fim)[2])
    for i in range(0, size(d_fim)[2]):
        for j in range(0, i):
            tigi[i,j] = trace_Matrix(ifim*d_fim[:,:,i]*ifim*d_fim[:,:,j])
            tigi[j,i] = tigi[i,j]
        
    d2 = d2fimdalpha2(alpha, model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpopdescr, ddescr, covd, sigma, docc, poped_db, 1e-4)["hess"].reshape(fim.size, hessp.size)
    d2logdfim = np.matmul(np.transpose(d2), np.asarray(ifim)).reshape(size(hessp)[0] ,size(hessp)[1])
    hess = -(hessp + d2logdfim - tigi)
    #   try({
    #     L=chol(fim)
    #     # calc inverse
    #     iL=solve(L,diag_matlab(length(L)))
    #     ifim=t(iL)%*%iL
    #     # calc trace of iF*dF/dAlpha(i)*iF*dF/dAlpha(j)
    #     for(i in 1:size(d_fim,3)){
    #       for(j in 1:i){
    #         tigi[i,j]=trace_Matrix(ifim*d_fim[,,i]*ifim*d_fim[,,j])
    #         tigi[j,i]=tigi[i,j]
    #       }
    #     }
    #     d2=d2fimdalpha2(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped_db,1e-4)
    #     d2Fim=reshape_matlab(d2,length(fim)^2,length(hessp)^2)
    #     d2logdfim=t(d2Fim)%*%ifim
    #     hess=-(hessp+reshape_matlab(d2logdfim,length(hessp),length(hessp))-tigi)
    #   })
    #   catch
    #   if((Engine$Type==1)){
    #     exception = lasterror
    #     if(exception$identifier=='MATLAB:posdef'){
    #       hess=zeros(length(alpha))
    #     } else {
    #       rethrow(exception)
    #     }
    #   } else {
    #     exception = lasterr
    #     if(exception$identifier==''){
    #       hess=zeros(length(alpha))
    #     } else {
    #       stop(sprintf(exception))
    #     }
    #   }
    
    return hess

