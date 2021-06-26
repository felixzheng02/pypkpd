"""


Author: Caiya Zhang, Yuchen Zheng
"""


from project.feval import feval
from project.zeros import zeros

#Hessian over eta, evaluated at eta_hat => laplace approximation possible
def hessian_eta_complex(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped_db,return_gradient=False):
  
    bAutomatic = False
    epsi0 = zeros(1, len(poped_db["parameters"]["notfixed_sigma"]))
    n = len(b_ind)
    hess= zeros(n)   # Memory for the Hessian matrix
    g = zeros(n,1)    # Memory for the gradient vector
    
    if bAutomatic:
        raise Exception("Automatic differentiation not yet implemented in Python version of pypkpd")
        #     if((poped_db["settings"]Engine$Type==2) ){#FreeMat
        #         stop(sprintf('Automatic differentiation is not available in PopED with FreeMat'))
        #     }
        #     b_init = hessianinit(b_ind)
        #     fg_init=feval(poped_db["model"]fg_pointer,x,a,bpop,b_init,bocc_ind)
        #      returnArgs =  feval(poped_db["model"]ferror_pointer,model_switch,xt_ind,fg_init,epsi0,poped_db) 
        # val = returnArgs[[1]]
        # poped_db = returnArgs[[1]]
        #     hess = val$hx
    else:
        h = 1E-04
        h2 = h * h
        for k in range(n):
            eta_plus = b_ind
            eta_plus[k] = eta_plus[k] + h * 1j
        
        if return_gradient:
            g_plus = feval(poped_db["model"]["fg_pointer"],x,a,bpop,eta_plus,bocc_ind)
            ff_plus = feval(poped_db["model"]["ferror_pointer"],model_switch,xt_ind,g_plus,epsi0,poped_db)
 ##没写！！！！ pixel image
            g[k] = Im(ff_plus)/h             # the kth gradient
        
        for l in range(k, n):                                    # Hessian (off-diagonal)
            eta_plus2 = eta_plus
            eta_plus2[l] = eta_plus2[l]+h
            g_plus = feval(poped_db["model"]["fg_pointer"],x,a,bpop,eta_plus2,bocc_ind)
            ff_plus = feval(poped_db["model"]["ferror_pointer"],model_switch,xt_ind,g_plus,epsi0,poped_db)
            
            eta_plus2[l]=eta_plus[l]-h
            g_plus = feval(poped_db["model"]["fg_pointer"],x,a,bpop,eta_plus2,bocc_ind)
            ff_minus = feval(poped_db["model"]["ferror_pointer"],model_switch,xt_ind,g_plus,epsi0,poped_db)
            
 ##没写！！！！ pixel image
            hess[k,l]=sum(Im(ff_plus-ff_minus)/h2/2)    # Hessian (central + complex step)
            hess[l,k]=hess[k,l]                           #Make hessian symmetric
        
    return [hess, g]
