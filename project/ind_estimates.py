"""


AUTHOR: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
import scipy as sp
from project.zeros import zeros
from project.diag_matlab import diag_matlab
from project.LinMatrixH import LinMatrixH
from project.ind_likelihood import ind_likelihood


#This is the function that should be minimized, w$r.t eta
def min_function(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped_db,c1,c2,c3,lC,det_res_var,b_ind,return_deriv=False):
    bAnalytic = False
    li = ind_likelihood(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped_db,lC,det_res_var,b_ind) #Individual log likelihood
    ret = c1 + c2 + 1/2 * np.transpose(b_ind) * c3 * b_ind - li
    if return_deriv == True:
        if bAnalytic == False:
            dret = zeros(d.size,1)
            heta = 1e-8
            for i in range(d.size):
                bind_plus = b_ind
                bind_plus[i] = b_ind(i)+heta
                li_plus = ind_likelihood(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped_db,lC,det_res_var,bind_plus) #Individual log likelihood
                dret[i] = (li_plus-li)/heta
        #}# else {
        #  dret = sum(dmodel_deta(tdata,cdata,theta,eta))
        #dret=c3*eta-dret
    return [ret, dret]



def opt_fun(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped_db,c1,c2,c3,lC,det_res_var,b_z):
    min_function(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped_db,c1,c2,c3,lC,det_res_var,b_z)

#Get the emperical bayes estimates for one individual
#Written for PopED by JN
def ind_estimates (data,bpop,d,sigma,start_bind:np.ndarray,bInter,bUDDLike,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped_db):
  
    b_i = np.transpose(start_bind)
    c1 = ((np.transpose(start_bind))/2*np.log(2*np.pi)).size
    c2 = 1/2*np.log(np.linalg.det(d))
    c3 = np.linalg.inv(d)
  
    if bInter == False and bUDDLike == False :#Calculate only one variance for all samples
        #eps = zeros(size(tdata,1),size(sigma,2))
        h = LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped_db)
        res_var = diag_matlab(diag_matlab(np.matmul(np.matmul(h,sigma),h)))
        lC = np.linalg.solve(np.linalg.cholesky(res_var),diag_matlab(res_var.size)) #Calculate cholesky factorization
        det_res_var = np.linalg.det(res_var)
    else:
        lC = 0
        det_res_var = 0
  
  
    #This can be changed to any optimization method in MATLAB or user written
    #Could scale b_i by b_i/sqrt(OM_i)*0.01 because fminsearch uses
    #max(5%*b_i,1E-05)
    #b_z = fminsearch(@(b_z) min_function(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped_db,c1,c2,c3,lC,det_res_var,b_z),b_i,optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',1000))
    b_z = sp.optim(b_i,opt_fun)
    opt_fun(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped_db,c1,c2,c3,lC,det_res_var,b_z)
    #b_z = optim(b_i,opt_fun)
    eta = b_z
  
    return eta 
