"""
## Evaluate the Fisher Information Matrix (FIM)
## 
## Compute the FIM given the model, parameters, design and methods defined in the 
## PopED database. Some of the arguments coming from the PopED database can be overwritten;  
## by default these arguments are \code{None} in the 
## function, if they are supplied then they are used instead of the arguments from the PopED database.
## 
## @param pypkpd_db A PopED database.
## @param fim.calc.type The method used for calculating the FIM. Potential values:
## \itemize{
## \item 0 = Full FIM.  No assumption that fixed and random effects are uncorrelated.  
## \item 1 = Reduced FIM. Assume that there is no correlation in the FIM between the fixed and random effects, and set these elements in 
## the FIM to zero. 
## \item 2 = weighted models (placeholder).
## \item 3 = Not currently used.
## \item 4 = Reduced FIM and computing all derivatives with respect to the standard deviation of the residual unexplained variation (sqrt(SIGMA) in NONMEM). 
## This matches what is done in PFIM, and assumes that the standard deviation of the residual unexplained variation is the estimated parameter
## (NOTE: NONMEM estimates the variance of the residual unexplained variation by default). 
## \item 5 = Full FIM parameterized with A,B,C matrices & derivative of variance. 
## \item 6 = Calculate one model switch at a time, good for large matrices. 
## \item 7 = Reduced FIM parameterized with A,B,C matrices & derivative of variance.
## }
## @param approx.method Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI
## @param FOCE.num Number individuals in each step of FOCE approximation method 
## @param bpop_val The fixed effects parameter values.  Supplied as a vector.
## @param d_full A between subject variability Matrix (OMEGA in NONMEM).
## @param docc_full A between occasion variability Matrix.
## @param sigma_full A residual unexplained variability Matrix (SIGMA in NONMEM).
## @param model_switch A Matrix that is the same size as xt, specifying which model each sample belongs to.
## @param ni A vector of the number of samples in each group.
## @param xt A Matrix of sample times.  Each row is a vector of sample times for a group.
## @param x A Matrix for the discrete design variables.  Each row is a group.
## @param a A Matrix of covariates.  Each row is a group.
## @param groupsize A vector of the number of individuals in each group.
## @param deriv_type A number indicating the type of derivative to use:
## \itemize{
## \item 0=Complex difference 
## \item 1=Central difference 
## \item 20=Analytic derivative (placeholder) 
## \item 30=Automatic differentiation (placeholder)
## }
## @param ... Other arguments passed to the function.
# @inheritParams Doptim
# @inheritParams create_poped_database
## 
## @return The FIM.
## 
## @family FIM
## @family evaluate_design
## @family evaluate_FIM
## 
## @export


## Author: Caiya Zhang, Yuchen Zheng
"""


from project.mftot import mftot
from project.param_set import param_set

def evaluate_fim(pypkpd_db,
                fim_calc_type=None,
                approx_method=None, 
                FOCE_num = None,
                bpop_val=None,
                d_full=None,
                docc_full=None,
                sigma_full=None,
                model_switch=None,
                ni=None,
                xt=None,
                x=None,
                a=None,
                groupsize=None,
                deriv_type = None,
                **kwargs):
  
  
    bpop_val = param_set(bpop_val, pypkpd_db, "parameters", "param_pt_val", "bpop")
    d_full = param_set(d_full, pypkpd_db, "parameters", "param_pt_val", "d")
    docc_full = param_set(docc_full, pypkpd_db, "parameters", "param_pt_val", "docc")
    sigma_full = param_set(sigma_full, pypkpd_db, "parameters", "param_pt_val", "sigma")
    
    #   if(is.null(model_switch)) model_switch = pypkpd_db$downsized.design$model_switch
    #   if(is.null(ni)) ni = pypkpd_db$downsized.design$ni
    #   if(is.null(xt)) xt = pypkpd_db$downsized.design$xt
    #   if(is.null(x)) x = pypkpd_db$downsized.design$x
    #   if(is.null(a)) a = pypkpd_db$downsized.design$a
    #   if(is.null(groupsize)) groupsize = pypkpd_db$downsized.design$groupsize
    #   
    model_switch = param_set(model_switch, pypkpd_db, "design", "model_switch")
    ni = param_set(ni, pypkpd_db, "design", "ni")
    xt = param_set(xt, pypkpd_db, "design", "xt")
    x = param_set(x, pypkpd_db, "design", "x")
    a = param_set(a, pypkpd_db, "design", "a")
    groupsize = param_set(groupsize, pypkpd_db, "design", "groupsize")
    
    if fim_calc_type is not None: 
        pypkpd_db["settings"]["iFIMCalculationType"] = fim_calc_type
    if approx_method is not None: 
        pypkpd_db["settings"]["iApproximationMethod"] = approx_method
    if FOCE_num is not None: 
        pypkpd_db["settings"]["iFOCENumInd"] = FOCE_num
    
    if deriv_type is not None: 
        pypkpd_db["settings"]["m1_switch"] = deriv_type
        pypkpd_db["settings"]["m2_switch"] = deriv_type
        pypkpd_db["settings"]["hle_switch"] = deriv_type
        pypkpd_db["settings"]["gradff_switch"] = deriv_type
        pypkpd_db["settings"]["gradfg_switch"] = deriv_type
    

    output = mftot(model_switch,groupsize,ni,xt,x,a,bpop_val,d_full,sigma_full,docc_full,pypkpd_db,kwargs)
    FIM = output["ret"]
    
    return FIM
