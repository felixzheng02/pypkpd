
"""
Function translated automatically using 'matlab.to.r()'

Author: Caiya Zhang, Yuchen Zheng
"""


import project.grad_all as grad_all
import project.zeros as zeros
import project.ind_estimates as ind_estimates

def grad_bpop (func,select_par,nout,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db,subset=poped_db["parameters"]["notfixed_bpop"], offdiag = False):
  #----------Model linearization with respect to pop parameters
  #
  # use helper function to check for/include EBEs
  #
    dx_dbpop = grad_all(func, select_par, nout, model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, d, sigma, docc, poped_db, subset=subset, noPopED = False, offdiag=offdiag)
  
    return dx_dbpop 

# helper for m2
def helper_v_EBE (model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db):
  
    if poped_db["settings"]["bCalculateEBE"]:
        #zeros(size(b_ind)[1],size(b_ind)[2])
        b_ind_x = ind_estimates(poped_db["mean_data"], bpop, d, sigma, t(b_ind),(poped_db["settings"]["iApproximationMethod"]==2), False, model_switch, xt_ind, x, a, b_ind, bocc_ind, poped_db)
    else:
        b_ind_x = b_ind
  
  
    vx = v(model_switch,xt_ind,x,a,bpop,b_ind_x,bocc_ind,d,sigma,docc,poped_db)[[1]]
    dim(vx) = c(prod(dim(vx)),1)
    return [vx, poped_db]


# helper for m1
def helper_LinMatrix (model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db):
  
    epsi0 = zeros(1,len(poped_db["parameters"]["notfixed_sigma"]))
  
    # create linearized model
    if (poped_db["settings"]["iApproximationMethod"]==0 or poped_db["settings"]["iApproximationMethod"]==3) : #FO, FOI
        b_ind=zeros(poped_db["parameters"]["NumRanEff"],1)
    
  
    if((poped_db["settings"]["bCalculateEBE"])):
        b_ind = ind_estimates(poped_db["mean_data"],bpop,d,sigma,t(b_ind),(poped_db["settings"]["iApproximationMethod"]==2),False,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped_db)
  
  
    g_p = feval(poped_db["model"]["fg_pointer"],x,a,bpop,b_ind,bocc_ind)
  
    returnArgs = feval(poped_db["model"]["ferror_pointer"],model_switch,xt_ind,g_p,epsi0,poped_db)
    ferror = returnArgs[[1]]
    
    if (poped_db["settings"]["iApproximationMethod"]==0 or poped_db["settings"]["iApproximationMethod"]==3 or (isempty(b_ind) and isempty(bocc_ind))) :
        #FO, FOI
        if (poped_db["settings"]["bUseSecondOrder"]):
            hess_eta = zeros(length(xt_ind),1)
            for o in 1:length(xt_ind):
                hessian_eta = hessian_eta_complex(model_switch[o],xt_ind[o],x,a,bpop,b_ind,bocc_ind,poped_db)
                hess_eta[o] = 1/2*trace_matrix(hessian_eta*d)
            ferror = ferror + hess_eta  
    else:
        #FOCE, FOCEI
        returnArgs = LinMatrixL(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped_db)
        l_plus = returnArgs[[1]]

        if((isempty(b_ind)) ){#No IIV present
            l_plus = zeros(size(xt_ind,1), 1)
        else:
            l_plus = l_plus%*%b_ind

    occ_add_plus = zeros(size(xt_ind,1), 1)
    if poped_db["parameters"]["NumOcc"] != 0:
        for m in 1:poped_db["parameters"]["NumOcc"]: 
            returnArgs <- LinMatrixL_occ(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,m,poped_db)
            l_plus_occ <- returnArgs[[1]]
            occ_add_plus = occ_add_plus + l_plus_occ*(bocc_ind[,m])
    ferror = ferror-(l_plus+occ_add_plus)
  
    return [ferror, poped_db]

