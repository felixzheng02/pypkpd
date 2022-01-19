"""
## Evaluate the expectation of determinant the Fisher Information Matrix (FIM)
## using the Laplace approximation.
## 
## Compute the expectation of the \code{det(FIM)} using the Laplace
## approximation to the expectation. Computations are made based on the model,
## parameters, distributions of parameter uncertainty, design and methods
## defined in the PopED database or as arguments to the function.
## 
## This computation follows the method outlined in Dodds et al, 
## "Robust Population Pharmacokinetic Experiment Design" JPP, 2005, equation 16.
## 
## Typically this function will not be run by the user.  Instead use \code{\link{evaluate.e.ofv.fim}}.
## 
## @param x The design parameters to compute the gradient on.
## @inheritParams evaluate.fim
## @inheritParams create.poped.database
## @inheritParams Doptim
## @param xtopto the sampling times
## @param xopto the discrete design variables
## @param optxt If sampling times are optimized
## @param opta If continuous design variables are optimized
## @param aopto the continuous design variables
## @param method If 0 then use an optimization routine translated from PopED code written in MATLAB to
##        optimize the parameters in the Laplace approximation.  If 1 then use \code{\link{optim}} to compute both
##        k and the hessian of k (see Dodds et al, JPP, 2005 for more information). If 2 then use \code{\link{fdHess}}
##        to compute the hessian.
## @param return_gradient Should the gradient be returned.
## @param ... Arguments passed through from other functions, does not pass anything to another function.
##   
## @return The FIM and the hessian of the FIM.
##   
## @family FIM
## @family E-family
## @export
## @keywords internal
# @importFrom nlme fdHess


## Author: Caiya Zhang, Yuchen Zheng
## right now function only works for normal and log-normal priors
"""

import numpy as np
import warnings

from numpy.core.fromnumeric import transpose
from project.size import size
from project.util import trans
from project.zeros import zeros
from project.mftot import mftot
from project.getfulld import getfulld
from project.graddetmf import graddetmf
from project.dfimdalpha import dfimdalpha
from project.diag_matlab import diag_matlab
from project.hesskalpha2 import hesskalpha2
from project.d2fimdalpha2 import d2fimdalpha2
from project.evaluate_fim import evaluate_fim
from project.trace_matrix import trace_matrix
from project.log_prior_pdf import log_prior_pdf
from project.line_search_uc import line_search_uc


def comp_grad_1(alpha, model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpopdescr, ddescr, covd, sigma, docc, poped_db, grad_p):
			returnArgs = dfimdalpha(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped_db,1e-6) 
			d_fim = returnArgs[0]
			fim = returnArgs[1]
			ifim = np.linalg.inv(fim)
			ifim.reshape(ifim.size, 1)
			d_fim.reshape(fim.size, len(grad_p))
			gradlogdfim = np.matmul(np.transpose(d_fim), ifim)
			grad_k = -(gradlogdfim + grad_p)

def calc_k(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,
                    ddescr,covd,sigma,docc,poped_db,Engine,return_gradient=False):
    bpop = bpopdescr[:,2]
    bpop[bpopdescr[:,1] != 0] = alpha[1:sum(bpopdescr[:,1]!=0)]
    d = ddescr[:,2]
    d[ddescr[:,1] == 4] = alpha[sum(bpopdescr[:,1]!=0):(alpha)]
    d = getfulld(d, covd)
    retargs = mftot(model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpop, d, sigma, docc, poped_db)
    fim = retargs["ret"]
    
    det_fim = np.linalg.det(fim)
    if det_fim < np.finfo(float).eps:
        warnings.warn("Determinant of the FIM is not positive")
        if return_gradient is True: 
            return {"k": np.nan, "grad_k": np.nan}
        return np.array({"k": np.nan})
    
    returnArgs = log_prior_pdf(alpha, bpopdescr, ddescr, return_gradient=return_gradient) 
    logp = returnArgs[0]
    if return_gradient is True:
        grad_p = returnArgs["grad"]
    k = -logp - np.log(det_fim)
    
    if return_gradient is True:
        
        # foo = tryCatch.W.E( comp_grad_1(alpha, model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpopdescr, ddescr, covd, sigma, docc, poped_db, grad_p) )
        # is.numeric(foo$value)
        # 
        # tryCatch.W.E( numDeriv::grad(function(x) calc_k(x,model_switch,groupsize,ni,xtoptn,xoptn,
        #                                                 aoptn,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
        #                                                 return_gradient=F),
        #                              alpha) )
        # 
        grad_k = comp_grad_1(alpha, model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpopdescr, ddescr, covd, sigma, docc, poped_db, grad_p)
        ## if not positive definite set grad_k=zeros(length(alpha),1)
        
        #tryCatch(log(det(fim)), warning = function(w) browser())
        return {"k": k, "grad_k": grad_k}
    
    
    return np.array({"k": k})


def ed_laplace_ofv(model_switch,groupsize,ni,xtopto,xopto,aopto,
                           bpopdescr,ddescr,covd,sigma,docc,poped_db,
                           method = 1,
                           return_gradient = False,
                           x = np.empty(),*argv):
    
    optxt=poped_db["settings"]["optsw"][1], 
    opta=poped_db["settings"]["optsw"][3],
  
    if any(ddescr[:,0] != 0 and ddescr[:,0] != 4):
        raise Exception("Only lognormal prior is supported for random effects!") 
    
    if any(bpopdescr[:,0] != 0 and 
            bpopdescr[:,0] !=4 and 
            #bpopdescr[:,0] != 2 and 
            bpopdescr[:,0] != 1):
        #stop(sprintf('Only uniform, normal and lognormal priors are supported for fixed effects!')) 
        raise Exception("Only normal and lognormal priors are supported for fixed effects!") 
        
    Engine = {"Type": 1, "Version": version["version_string"]}
    
    x2 = x
    if len(x) != 0:
        if optxt is True:
            notfixed=poped_db["design_space"]["minxt"] != poped_db["design_space"]["maxxt"]
            if poped_db["design_space"]["bUseGrouped_xt"]:
                xtopto[notfixed] = x[poped_db["design_space"]["G_xt"][notfixed]]
                u_tmp, idx = np.unique(poped_db["design_space"]["G_xt"][notfixed], return_index=True)
                u = u_tmp[idx.argsort()]
                x[1:np.prod(size(u))] = np.zeros(1)
            else:
                xtopto[notfixed]=x[0:np.prod(size(xtopto[notfixed]))]
                x = x[-np.arrange(1, np.prod(size(xtopto[notfixed])))]
        
        if opta is True:
            notfixed = poped_db["design_space"]["mina"] != poped_db["design_space"]["maxa"]
            if poped_db["design_space"]["bUseGrouped_a"] is True:
                aopto[notfixed] = x[poped_db["design_space"]["G_a"][notfixed]]
            else:
                aopto[notfixed] = x
        x = x2    
    
    
    # alpha parameter vector
    alpha_k = np.array([bpopdescr[bpopdescr[:,0]!=0, 2], ddescr[ddescr[:,0]!=0, 2]]).reshape(2,1)
    
    ## do log transformation of ln bpop and d parameters for unconstrained optimization 
    if alpha_k.size > sum(bpopdescr[:,0] != 0):
        d_index = range(sum(bpopdescr[:,0] != 0), len(alpha_k))
    else:
        d_index = None
    
    bpop_index = range(0, sum(bpopdescr[:,0] != 0))
    unfixed_bpop = bpopdescr[bpopdescr[:,0] != 0,:]
    exp_index = np.array([unfixed_bpop[:,0] == 4, d_index == d_index])
    alpha_k_log = alpha_k
    if any(exp_index) is True: 
        alpha_k_log[exp_index] = np.log(alpha_k[exp_index])
    #alpha_k_log[d_index]=log(alpha_k[d_index])


    
    #calc initial k value and gradient
    ret_grad = False
    if method == 0: 
        ret_grad = True
    
    ret_args = calc_k(alpha_k,model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
                        return_gradient=ret_grad) 
    f_k = ret_args[0]
    gf_k = None
    if method == 0:
        gf_k = ret_args["grad_k"]
    
    ret_args = tryCatch.W.E(np.linalg.inv(evaluate_fim(poped_db)))
    if any(type(ret_args["value"]) is "error"):
        warnings.warn("Inversion of the FIM is not possible, the current design cannot estimate all parameters")
        warnings.warn("In "+ str(ret_args["value"]["call"]) + " : \n    " + ret_args["value"]["message"])
        f = np.nan
        if return_gradient: 
            return {"f": f, "gf": np.nan}

        return {"f": f}
    # tryCatch.W.E( numDeriv::grad(function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
    #                                                 aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
    #                                                 return_gradient=F),
    #                              alpha_k) )
    # 
    # tryCatch.W.E( numDeriv::hessian(function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
    #                                                 aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
    #                                                 return_gradient=F),
    #                              alpha_k) )
    # 
    #   ## test f_k value
    #fim = det(evaluate.fim(poped_db))
    #   # assuming all normal distributions and only first 3 fixed effects
    #p = prod(dnorm(bpopdescr[1:3,2],mean=bpopdescr[1:3,2],sd=sqrt(bpopdescr[1:3,3]))) 
    #k_test = -log(fim*0.36)
    #   k_test==f_k
    
    #   ## test gf_k
    #   alpha_k_plus = alpha_k+rbind(0.00001,0,0)
    #   returnArgs.1 = calc_k(alpha_k_plus,model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
    #                        return_gradient=F) 
    #   f_k_plus = returnArgs.1[[1]]
    # 
    #   gf_k_test.1  = (f_k_plus - f_k)/0.00001
    
    #transform gradient for ds (log(alpha))
    #   if(!isempty(d_index)){
    #     gf_k[d_index]=gf_k[d_index]*exp(alpha_k_log[d_index])
    #   }
    if exp_index.size != 0 and gf_k is not None:
        gf_k[exp_index] = gf_k[exp_index]*np.exp(alpha_k_log[exp_index])
    
    if np.isnan(f_k):
        f = 0
        gf = zeros(size(x))
        if return_gradient is True: 
            return {"f": f, "gf": gf}

        return {"f": f}
    
    
    ###################
    ## minimization of k(alpha) 
    ###################
    
    if method == 0: ## sebastian method
        #initialize optimization variables
        dim = alpha_k.size
        H_k = diag_matlab(dim)
        B_k = H_k
        niter = 0
        while np.linalg.norm(gf_k, type="2") > 0.001: # while inner conv. krit. not met
            #determine search direction for line search
            p_k = -np.matmul(H_k, gf_k)
            f_name  = "calc_k"
            #f_options = list(trans(alpha),model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine)
            f_options = ["replace",model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine]
            returnArgs = line_search_uc(alpha_k_log,f_k,gf_k,p_k,f_name,f_options,exp_index)
            alpha_k1_log = returnArgs[0]
            f_k1 = returnArgs[1]
            s_k = alpha_k1_log - alpha_k_log
            if max(abs(np.transpose(s_k))/max(np.array([np.transpose(alpha_k1_log), np.ones(1).reshape(1,np.transpose(alpha_k1_log).size)]).reshape(2,1))) < 1e-3:  
                # check this that it is the same as in matlab
                break
            
            f_k1 = calc_k(trans(alpha_k1_log),model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,return_gradient=True) 
            gf_k1 = attr(f_k1,"grad")
            #transform gradient for ds (log(alpha))
            #     if(!isempty(d_index)){
            #       gf_k1[d_index]=gf_k1[d_index]*exp(alpha_k1_log(d_index))
            #     }
            if exp_index.size != 0:
                gf_k1[exp_index] = gf_k1[exp_index]*np.exp(alpha_k1_log[exp_index])
            
            y_k = gf_k1 - gf_k
            rho_k = 1/np.matmul(np.transpose(y_k), s_k)
            rho_k = rho_k[:,:]
            if np.matmul(np.transpose(y_k), s_k)/-np.matmul(np.transpose(gf_k), s_k) > np.finfo(float).eps:
                tmp1 = np.matmul(rho_k*s_k, np.transpose(y_k))
                tmp2 = np.matmul(rho_k*y_k, np.transpose(s_k))
                H_k= np.matmul(np.matmul((diag_matlab(dim) - tmp1), H_k), (diag_matlab(dim)-tmp2)) + np.matmul(rho_k*s_k, np.transpose(s_k))
      

            alpha_k_log = alpha_k1_log
            gf_k = gf_k1
            f_k = f_k1
            niter = niter + 1
        
        alpha_k = trans(alpha_k_log)
        
        #if the number of iterations is smaller than the dimension of the problem
        #we have to calculate the hessian explicitly
        if niter < B_k.size or poped_db["settings"]["iEDCalculationType"] == 1:
            hess = hesskalpha2(alpha_k, model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,1e-6,Engine)
            detHessPi = np.linalg.det(hess)*(2*np.pi)^(-hess.size)
        else:
            temp = np.ones(size(gf_k)).reshape(size(gf_k),1)
            #temp[d_index]=1/alpha_k[d_index]
            temp[exp_index] = 1/alpha_k[exp_index]
            
            iH_k = np.linalg.inv(H_k)
            hess = iH_k*np.matmul(temp, np.transpose(temp))
            detHessPi = np.linalg.det(hess)*(2*np.pi)^(-hess.size)
    else:  # end sebastian method
        
        #     priordescr = rbind(bpopdescr,ddescr)
        #     priordescr = priordescr[priordescr[,1]==2,]
        #     lb = priordescr[,2]-priordescr[,3]/2
        #     ub = priordescr[,2]+priordescr[,3]/2
        
        
        ## minimize K(alpha_k)
        opt_meth = "Nelder-Mead"
        if alpha_k.size == 1: 
            opt_meth = "BFGS"
        # initial_k = calc_k(alpha_k,model_switch,groupsize,ni,xtopto,xopto,
        #                     aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
        #                     return_gradient=F)
        # cat("alpha_k = ",alpha_k," initial k = ", initial_k, "\n")
        # cat("xt = ",xtopto, "\n")
        # 
        # bpop=bpopdescr[,2,drop=F]
        # bpop[bpopdescr[,1,drop=F]!=0]=alpha_k[1:sum(bpopdescr[,1,drop=F]!=0),drop=F]
        # d=ddescr[,2,drop=F]
        # d[ddescr[,1]==4]=alpha_k[(sum(bpopdescr[,1,drop=F]!=0)+1):length(alpha_k),drop=F]
        # d=getfulld(d,covd)
        # retargs=mftot(model_switch,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db)
        # fim = retargs$ret
        # det(fim)
        # inv(fim)
        # 
        # for(tmp in seq(0.0993,10, by=1)){
        #   cat(calc_k(tmp,model_switch,groupsize,ni,xtopto,xopto,
        #              aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
        #              return_gradient=F),"\n")
        # }
        # 
        # nlme::fdHess(alpha_k,
        #              function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
        #                                 aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
        #                                 return_gradient=F)) 
        # 
        # numDeriv::hessian(function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
        #                                      aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
        #                                      return_gradient=F),
        #                   alpha_k)
        # 
        output = optim(alpha_k, 
                        function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
                                        aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
                                        return_gradient=F),
                        #gr=function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
                        #                     aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
                        #                     return_gradient=T)[["grad_k"]],
                        method=opt_meth,
                        #method="BFGS",
                        #method="Brent",
                        #lower=1e-6,
                        #upper=1e6,
                        #lower=0.0000001,
                        #lower=lb,
                        #upper=1000,
                        hessian=T)
        hess = output["hessian"]
        f_k = output["value"]
        
        ## Different hessian and gradient calculation
        if method == 2: 
            if requireNamespace("nlme", quietly = True) is not True:
                raise Exception("nlme package needed for this function to work with option 'method=2'. Please install it.")
            
        k_vals = nlme::fdHess(output["par"],
                                function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
                                                    aopto,bpopdescr,ddescr,covd,sigma,docc,poped_db,Engine,
                                                    return_gradient=F)) 
        hess = k_vals["Hessian"]
        
    #f=Re(-exp(-f_k)/sqrt(detHessPi))
    det_hess_pi = np.linalg.det(hess*(2*np.pi)^(-1))
    if det_hess_pi < 0:
        warnings.warn("The laplace OFV is nan, because det(hessian) < 0.")
        f = np.nan
    else:
        f = np.sqrt(det_hess_pi)^(-1)*np.exp(-f_k)
    
    if return_gradient is True:
        bpop = bpopdescr[:,1]
        bpop[bpopdescr[:,0] != 0] = alpha_k[1:sum(bpopdescr[:,0] != 0)]
        d = ddescr[:,1]
        d[ddescr[:,0] != 0] = alpha_k[sum(bpopdescr[:,0] != 0):alpha_k.size]
        d = getfulld(d, covd)
        
        gradxt = np.zeros(1)
        grada = np.zeros(1)
        if optxt is True:
            notfixed = poped_db["design_space"]["minxt"] != poped_db["design_space"]["maxxt"]
            gradxt = -graddetmf(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True,gradxt=True)
            gradxt = gradxt(notfixed)
            if poped_db["design_space"]["bUseGrouped_xt"]:                
                index_tmp, idx = np.unique(poped_db["design_space"]["G_xt"], return_index=True) 
                index = index_tmp[idx.argsort()] 
                gradxt=gradxt(index)
            
        if opta is True:
            notfixed = poped_db["design_space"]["mina"]!=poped_db["design_space"]["maxa"]
            grada = -graddetmf(model_switch, np.ones(size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True)
            grada = grada(notfixed)
            if poped_db["design_space"]["bUseGrouped_a"]:
                index_tmp, idx = np.unique(poped_db["design_space"]["G_a"], return_index=True) 
                index = index_tmp[idx.argsort()]
                grada = grada(index)
            
        dkdxt = np.array([gradxt,grada]).reshape(1,2)
        
        h_alpha = 1e-4
        tensor = np.zeros(1).reshape(alpha_k.size, alpha_k.size, dkdxt.size)
        for i in range(0, alpha_k.size):
            for j in range(0, i):

                alpha_plus_plus = alpha_k
                alpha_plus_plus[i] = alpha_plus_plus[i] + h_alpha
                alpha_plus_plus[j] = alpha_plus_plus[j] + h_alpha
                bpop = bpopdescr[:,1]
                bpop[bpopdescr[:,0] != 0] = alpha_plus_plus[1:sum(bpopdescr[:,0] != 0)]
                d = ddescr[:,1]
                d[ddescr[:,0] != 0] = alpha_plus_plus[sum(bpopdescr[:,0] != 0):alpha_plus_plus.size]
                d = getfulld(d, covd)
                if optxt is True:
                    notfixed = poped_db["design_space"]["minxt"] != poped_db["design_space"]["maxxt"]
                    gradxt = np.transpose(graddetmf(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True,gradxt=True))
                    gradxt = gradxt(notfixed)
                    if poped_db["design_space"]["bUseGrouped_xt"]:
                        index_tmp, idx = np.unique(poped_db["design_space"]["G_xt"], return_index=True) 
                        index = index_tmp[idx.argsort()]
                        gradxt = gradxt(index)
                    
                if opta is True:
                    notfixed = poped_db["design_space"]["mina"] != poped_db["design_space"]["maxa"]
                    grada = -graddetmf(model_switch,np.ones(size(aopto)).reshape(size(aopto),1),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True)
                    grada = grada(notfixed)
                    if poped_db["design_space"]["bUseGrouped_a"]:
                        index_tmp, idx = np.unique(poped_db["design_space"]["G_a"], return_index=True) 
                        index = index_tmp[idx.argsort()]
                        grada = grada(index)
                    
                dkdxt_plus_plus = np.array([gradxt,grada])
                
                alpha_minus_plus = alpha_k
                alpha_minus_plus[i] = alpha_minus_plus[i] - h_alpha
                alpha_minus_plus[j] = alpha_minus_plus[j] + h_alpha
                bpop = bpopdescr[:,1]
                bpop[bpopdescr[:,0] != 0] = alpha_minus_plus[0:sum(bpopdescr[:,0] != 0)]
                d = ddescr[:,1]
                d[ddescr[:,0] != 0] = alpha_minus_plus[sum(bpopdescr[:,0] != 0):alpha_minus_plus.size]
                d = getfulld(d, covd)
                if optxt is True:
                    notfixed = poped_db["design_space"]["minxt"] != poped_db["design_space"]["maxxt"]
                    gradxt = np.transpose(graddetmf(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True,gradxt=True))
                    gradxt = gradxt(notfixed)
                    if poped_db["design_space"]["bUseGrouped_xt"]:
                        index_tmp, idx = np.unique(poped_db["design_space"]["G_xt"], return_index=True) 
                        index = index_tmp[idx.argsort()]
                        gradxt = gradxt(index)
                
                if opta is True:
                    notfixed = poped_db["design_space"]["mina"] != poped_db["design_space"]["maxa"]
                    grada = -graddetmf(model_switch,np.ones(size(aopto)).reshape(size(aopto),1),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True)
                    grada = grada(notfixed)
                    if poped_db["design_space"]["bUseGrouped_a"]:
                        index_tmp, idx = np.unique(poped_db["design_space"]["G_a"], return_index=True) 
                        index = index_tmp[idx.argsort()]
                        grada = grada(index)
                    
                dkdxt_minus_plus = np.array([gradxt,grada])
                
                alpha_plus_minus = alpha_k
                alpha_plus_minus[i] = alpha_plus_minus[i] + h_alpha
                alpha_plus_minus[j] = alpha_plus_minus[j] - h_alpha
                bpop = bpopdescr[:,1]
                bpop[bpopdescr[:,0] != 0] = alpha_minus_plus[1:sum(bpopdescr[:,0] != 0)]
                d = ddescr[:,1]
                d[ddescr[:,0] != 0] = alpha_minus_plus[sum(bpopdescr[:,0] != 0):alpha_minus_plus.size]
                d = getfulld(d, covd)
                if optxt is True:
                    notfixed = poped_db["design_space"]["minxt"] != poped_db["design_space"]["maxxt"]
                    gradxt = np.transpose(graddetmf(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True,gradxt=True))
                    gradxt = gradxt(notfixed)
                    if poped_db["design_space"]["bUseGrouped_xt"]:
                        index_tmp, idx = np.unique(poped_db["design_space"]["G_xt"], return_index=True) 
                        index = index_tmp[idx.argsort()]
                        gradxt = gradxt(index)
                
                if opta is True:
                    notfixed = poped_db["design_space"]["mina"] != poped_db["design_space"]["maxa"]
                    grada = -graddetmf(model_switch,np.ones(size(aopto)).reshape(size(aopto),1),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True)
                    grada = grada(notfixed)
                    if poped_db["design_space"]["bUseGrouped_a"]:
                        index_tmp, idx = np.unique(poped_db["design_space"]["G_a"], return_index=True) 
                        index = index_tmp[idx.argsort()]
                        grada = grada(index)
                
                dkdxt_plus_minus = np.array([gradxt,grada])
                
                alpha_minus_minus = alpha_k
                alpha_minus_minus[i] = alpha_minus_minus[i] - h_alpha
                alpha_minus_minus[j] = alpha_minus_minus[j] - h_alpha
                bpop = bpopdescr[:,1]
                bpop[bpopdescr[:,0] != 0] = alpha_minus_plus[1:sum(bpopdescr[:,0] != 0)]
                d = ddescr[:,1]
                d[ddescr[:,0] != 0] = alpha_minus_plus[sum(bpopdescr[:,0] != 0):alpha_minus_plus.size]
                d = getfulld(d, covd)
                if optxt is True:
                    notfixed = poped_db["design_space"]["minxt"] != poped_db["design_space"]["maxxt"]
                    gradxt = np.transpose(graddetmf(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True,gradxt=True))
                    gradxt = gradxt(notfixed)
                    if poped_db["design_space"]["bUseGrouped_xt"]:
                        index_tmp, idx = np.unique(poped_db["design_space"]["G_xt"], return_index=True) 
                        index = index_tmp[idx.argsort()]
                        gradxt = gradxt(index)
                
                if opta is True:
                    notfixed = poped_db["design_space"]["mina"] != poped_db["design_space"]["maxa"]
                    grada = -graddetmf(model_switch,np.ones(size(aopto)).reshape(size(aopto),1),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped_db,lndet=True)
                    grada = grada(notfixed)
                    if poped_db["design_space"]["bUseGrouped_a"]:
                        index_tmp, idx = np.unique(poped_db["design_space"]["G_a"], return_index=True) 
                        index = index_tmp[idx.argsort()]
                        grada = grada(index)
                    
                dkdxt_minus_minus = np.array([gradxt,grada])
            
            tensor[i,j,:] = ((dkdxt_plus_plus-dkdxt_plus_minus-dkdxt_minus_plus+dkdxt_minus_minus))/(4*h_alpha^2)
            tensor[j,i,:] = tensor[i,j,:]
        
        ddetHessdxt = zeros(dkdxt.size, 1)
        for i in range(0, dkdxt.size):
            ddetHessdxt[i] = detHessPi*trace_matrix(np.linalg.inv(hess)*(-tensor[:,:,i]))
        
        
        gf = Re(np.exp(-f_k)*(2*detHessPi*dkdxt+ddetHessdxt)/(2*detHessPi^(3/2)))
        
        return {"f": f, "gf": gf}
    
    return {"f": f}
    




