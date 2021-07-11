"""
## Function translated automatically using 'matlab.to.r()'


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from project.feval import do_call
from project.grad_bpop import grad_bpop


def m3(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,bUseVarSigmaDerivative,poped_db,helper_v_EBE):
    #
    # size: (samps per subject^2 x (number of random effects + number of occasion variances + number of sigmas))

    dv_dd = None
    dv_covd = None
    dv_ddocc = None
    dv_covdocc = None
    dv_dsig = None
    dv_dcovsig = None

    if sum(poped_db["parameters"]["notfixed_d"]) > 0:
        dv_dd = grad_bpop(helper_v_EBE,8,size(xt_ind)[0]^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db, subset=poped_db["parameters"]["notfixed_d"])
  

    if sum(poped_db["parameters"]["notfixed_covd"]) > 0:
        dv_covd = grad_bpop(helper_v_EBE,8,size(xt_ind)[0]^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db, subset=poped_db["parameters"]["notfixed_covd"], offdiag = True)
    

    if sum(poped_db["parameters"]["notfixed_docc"]) > 0:
        dv_ddocc = grad_bpop(helper_v_EBE,10,size(xt_ind)[0]^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db, subset=poped_db["parameters"]["notfixed_docc"])
 

    if sum(poped_db["parameters"]["notfixed_covdocc"]) > 0:
        dv_covdocc = grad_bpop(helper_v_EBE,10,size(xt_ind)[0]^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db, subset=poped_db["parameters"]["notfixed_covdocc"], offdiag = True)
  

    if sum(poped_db["parameters"]["notfixed_sigma"]) > 0:
        dv_dsig = grad_bpop(helper_v_EBE,9,size(xt_ind)[0]^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db, subset=poped_db["parameters"]["notfixed_sigma"])
        if bUseVarSigmaDerivative == False:
            dv_dsig = np.transpose(2 * np.sqrt(np.diag(sigma))[poped_db["parameters"]["notfixed_sigma"]==1] * np.transpose(dv_dsig))


    if sum(poped_db["parameters"]["notfixed_covsigma"]) > 0:
        dv_dcovsig = grad_bpop(helper_v_EBE,9,size(xt_ind)[0]^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db, subset=poped_db["parameters"]["notfixed_covsigma"], offdiag = True)
  

    dv_db = do_call(np.hstack, [dv_dd, dv_covd, dv_ddocc, dv_covdocc, dv_dsig, dv_dcovsig])

    return {"dv_db": dv_db, "poped_db": poped_db}
