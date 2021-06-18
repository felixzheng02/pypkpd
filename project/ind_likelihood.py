"""#Calculates the individual ln likelihood
# function li=posthoc_dist(data,xt,x,a,bpop,d,bhat,...
#                          EPS,model_switch,...
#                          hlf,hlg,hle,...
#                          INTER)
# % from davidian and giltinan p 173 
# dmat= diag(d);
# 
# fg_bhat=fg(x,a,bpop,bhat);
# ff_bhat = ff(model_switch,xt,fg_bhat);
# 
# res = data - ff_bhat;
# 
# if INTER==true
# h=LinMatrixH_foce(model_switch,xt,x,a,bpop,EPS,hle,bhat);
# else
#   h=LinMatrixH(model_switch,xt,x,a,bpop,EPS,hle);
# end
# RR = diag(diag(h*diag(EPS)*h'));
# 
#   %li = log(det(dmat))+ (bhat'/dmat)*bhat+ ...
#           li = (bhat'/dmat)*bhat+ ...
#        log(det(RR))+ (res'/RR)*res;        


Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.feval import feval
from project.LinMatrixH import LinMatrixH
from project.diag_matlab import diag_matlab


def ind_likelihood(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped_db,lC,det_res_var,b_ind):
    # % from davidian and giltinan p 173 
    if bUDDLike == False:
        #browser()
        fg_bhat = feval(poped_db["model"]["fg_pointer"],x,a,bpop,b_ind,bocc_ind)
        ipred = feval(poped_db["model"]["ff_pointer"],model_switch,xt_ind,fg_bhat,poped_db)[["y"]]
        res = data-ipred#Individual residuals
    
        if bInter == True:
            #For Cases WITH interaction, linearize around eta = eta^
            #eps = zeros(size(tdata),1),size(sigma,2))
            #eps = zeros(size(t(data),1),size(sigma,2))
            h = LinMatrixH(np.transpose(model_switch),
                     np.transpose(xt_ind),
                     np.transpose(x),
                     np.transpose(a),
                     bpop,
                     b_ind,
                     bocc_ind,
                     poped_db)["y"] #The covariance for this individual
            res_var = diag_matlab(diag_matlab(np.transpose(np.matmul(h,sigma),np.transpose(h))))
            lC = np.linalg.solve(np.linalg.cholesky(res_var),diag_matlab(res_var.size))
            det_res_var = np.linalg.det(res_var)
        #else:
            #Cases WITHOUT interaction, linearize around eta = 0
            #  h = LinMatrixH(tdata,cdata,theta,zeros(size(eta)),eps) #The covariance for this individual
    
        R = (np.np.matmul(np.transpose(res),lC))
        li = -1/2*np.log(det_res_var)-1/2*np.matmul(R,np.transpose(R)) # + const
    
    else:
        #%UDD likelihood
        #li=sum(model(tdata,cdata,theta,eta))
        raise Exception("User defined likelihood not implemented for PopED in R")  
  
    return li

