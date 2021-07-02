"""



Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.size import size
from project.zeros import zeros
from project.diag_matlab import diag_matlab
from project.get_all_params import get_all_params
from project.get_unfixed_params import get_unfixed_params

def get_cv(param_vars,poped_db):
    #Return the RSE,CV of parameters
    ## Author: Andrew Hooker
    params_all =  get_all_params(poped_db)[[7]] 
    
    returnArgs =  get_unfixed_params(poped_db,params_all) 
    params = returnArgs[[7]]
    var_derivative = returnArgs[[8]]
    
    if param_vars.size != params.size:
        raise Exception("Number of unfixed parameters not the same as the size of the FIM,\nno RSE can be computed!\n")
    
    params_cv = zeros(size(param_vars))
    for i in range(0, params.size):
        if params[i] != 0:
            if var_derivative[i] == 1:
                params_cv[i] = np.sqrt(param_vars[i])/abs(params[i])
            else:   #Derivative w$r.t to SD instead of var
                params_cv[i] = np.sqrt(param_vars[i])/np.sqrt(params[i])
        
        else:
            params_cv[i] = np.sqrt(param_vars[i])
    
    return [params, params_cv]


#' Compute the expected parameter relative standard errors 
#' 
#' This function  computes the expected relative standard errors of a model given a design and a previously computed
#' FIM.
#' 
#' @param fim A Fisher Information Matrix (FIM).
#' @param bpop A vector containing the values of the fixed effects used to compute the \code{fim}. 
#' @param d A vector containing the values of the diagonals of the between subject variability matrix.
#' @param use_percent Should RSE be reported as percent? 
#' @param prior_fim A prior FIM to be added to the \code{fim}. Should be the same size as the \code{fim}.
#' @param ... Additional arguments passed to \code{\link{inv}}. 
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' 
#' @return A named list of RSE values.  If the estimated parameter is assumed to be zero then for that 
#'   parameter the standard error is returned.
#' 
#' @family evaluate_design
#' 
# @example inst/examples_fcn_doc/examples_evaluate.fim.R
#' @example tests/testthat/examples_fcn_doc/examples_evaluate.fim.R
#' @export

def get_rse(fim, poped_db,*args):
    bpop = poped_db["parameters"]["bpop"][:,1]
    #bpop=poped_db["parameters"]bpop[,2,drop=F],
    d = poped_db["parameters"]["d"][:,1]
    # d = poped_db["parameters"]d[,2,drop=F],
    docc = poped_db["parameters"]["docc"]
    sigma = poped_db["parameters"]["sigma"]
    use_percent = True
    fim_calc_type = poped_db["settings"]["iFIMCalculationType"],
    prior_fim = poped_db["settings"]["prior_fim"],
    #pseudo_on_fail = False,
    ## update poped_db with options supplied in function
    called_args = match_call()
    default_args = formals()
    for i in called_args.keys()[-1]:
        if length(grep("^poped\\.db\\$",capture_output(default_args[[i]]))) == 1:
        #eval(parse(text=paste(capture.output(default_args[[i]]),"=",called_args[[i]])))
        # if (i %in% c('bpop','d')) {
        #   if (eval(parse(text=paste("dim(",i,")[2]>1"))))
        #     (eval(parse(text=paste(i, "=",i,"[,2]"))))
        # }
            eval(parse(text=paste(capture_output(default_args[[i]]),"=",i)))

    ## if prior is given in poped_db then add it to the given fim
    if len(prior_fim) != 0 and all(size(prior_fim) == size(fim)):
        fim = fim + prior_fim
  
  
    inv_fim = tryCatch({inv(fim,...)}, error=function(e){
        warning(e)
        return None
    })

    if inv_fim is None:
        mess = paste0("\n  Could not invert the FIM.",
                    "\n  Is the design adequate to estimate all parameters?")
        eig = eigen(fim)[["values"]]
        names(eig) = get_parnam(poped_db)
        neg_vals = eig[eig< 0]
        num_neg = length(neg.vals)
        if num_neg > 0:
            mess = paste0(mess,"\n  Potentially problematic parameters and associated eigenvalues:")
            for i in range(0,num_neg):
                mess = paste0(mess,sprintf("\n %12s  %8.7e",names(neg.vals[i]),neg_vals[i]))
        #warning(simpleWarning(mess,call="get_rse()"))
        warning(mess)
        return (np.repeat(np.nan, length(get_parnam(poped_db))))

    param_vars = diag_matlab(inv_fim)
    returnArgs =  get_cv(param_vars,poped_db) 
    params = returnArgs[[1]]
    params_rse = returnArgs[[2]]
    parnam = get_parnam(poped_db)
    ret = params_rse[:,:,drop=T]
    if use_percent: 
        ret[params!=0] = ret[params!=0]*100
    ret.keys() = parnam
    if any(ret==0):
        zero_ret = ret[ret==0].keys()
        mess = paste0("  The following parameters are not estimable:\n  ",
                    paste0(zero_ret,collapse = ", "),
                    "\n  Is the design adequate to estimate all parameters?")
        warning(mess, call. = False)
        ret[ret==0] = np.nan
    
    return ret

