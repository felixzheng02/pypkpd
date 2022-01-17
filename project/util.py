"""
is_not_none(d: dict, k: str) -> boolean
Judges if a str is in a dict's key and has corresponding value.
Returns True if it is in the dict.
Returns False if it is not.
@param d: dictionary {key, value}
@param k: the string needs to be searched in the key.
@return: boolean

dots

unbound_par

bound_par

transform_back_par

get_fim_size(dict) -> int
Makes some calculations and returns an int based on a poped database input.
@param poped_db: poped database as a dictionary
@return: "fim_size" as an integer based on some calculations

Author: Caiya Zhang, Yuchen Zheng
"""

import numpy as np
import scipy as sp

def is_not_none(d: dict, k: str):
	if k in list(d.keys()):
		if d[k] is not None:
			return True
		else:
			return False
	else:
		return False

def trans(x, exp_index):
        x[exp_index] = np.exp(x[exp_index])
        return x
#trans  = function(x) matrix(c(x[bpop_index],exp(x[d_index])),ncol=1,byrow=T


def dots(*args):
  eval(substitute(alist(*args)))


def unbound_par(x,lower=-np.Infinity, upper=np.Infinity,tol=1e-8,logit=True):
  
	if upper < lower:
		raise Exception("'upper' must be greater than 'lower.'")
	if all(np.isfinite(np.array([lower,upper]))) and (tol > (upper - lower)/2):
		raise Exception("'tol' must be less than half the distance between the upper and lower bounds.")
	
	lower_tol:np.ndarray = lower + tol
	upper_tol:np.ndarray = upper - tol
	x_tmp = np.stack([x, lower_tol])
	x = x_tmp.max(0)
	x_tmp = np.stack([x, upper_tol])
	x = x_tmp.min(0)
	
	if all(np.isinf(np.array([lower,upper]))):
		return x
	elif np.isinf(lower):
		return np.log(upper - x)
	elif np.isinf(upper):
		return np.log(x - lower)
	elif logit is True:
		return sp.stats.norm.logppf((x - lower)/(upper-lower))
	else: # probit
		return sp.stats.norm.ppf((x - lower)/(upper-lower))
	

def bound_par(x,lower=-np.Infinity, upper=np.Infinity,logit=True):
  
	if upper < lower: 
		raise Exception("'upper' must be greater than 'lower.'")
	
	if all(np.isinf(np.array([lower,upper]))):
		return(x)
	elif np.isinf(lower):
		return upper - np.exp(x)
	elif np.isinf(upper):
		return np.exp(x) + lower
	elif logit is True:
		return sp.stats.norm.logpdf(x)*(upper - lower) + lower
	else: # probit
		return sp.stats.norm.pdf(x)*(upper - lower) + lower
	


# transform_back = function(par,lower=-Inf,upper=Inf){
#   # FastImputation::BoundNormalizedVariable(
#   #   par,
#   #   constraints = 
#   #     list(lower=lower,
#   #          upper=upper))
#   bound_par(par,lower=lower,upper=upper)
# }

def transform_back_par(ps_tbl,*args):

	par=ps_tbl["par"],
	ps_transformed=ps_tbl["transformed"],
	ps_lower_orig=ps_tbl["lower_orig"],
	ps_upper_orig=ps_tbl["upper_orig"]

	if any(ps_transformed) is True:
		par[ps_transformed] = mapply(bound_par, par[ps_transformed], ps_lower_orig[ps_transformed], ps_upper_orig[ps_transformed])
	
	return par 


"""
##' Catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
##'
##' @title tryCatch both warnings (with value) and errors
##' @param expr an \R expression to evaluate
##' @return a list with 'value' and 'warning', where
##'   'value' may be an error caught.
##' @author Martin Maechler, The R Core Team
##' @keywords internal
def tryCatch_W_E(expr):
	W = None
	w_handler = function(w){ # warning handler
		W = w
		invokeRestart("muffleWarning")
	}
	return {"value": withCallingHandlers(tryCatch(expr, error = function(e) e),warning = w.handler),
			"warning": W}
"""


def get_fim_size(poped_db):
	numnotfixed_bpop = sum(poped_db["parameters"]["notfixed_bpop"])
	numnotfixed_d    = sum(poped_db["parameters"]["notfixed_d"])
	numnotfixed_covd = sum(poped_db["parameters"]["notfixed_covd"])
	numnotfixed_docc  = sum(poped_db["parameters"]["notfixed_docc"])
	numnotfixed_covdocc  = sum(poped_db["parameters"]["notfixed_covdocc"])
	numnotfixed_sigma  = sum(poped_db["parameters"]["notfixed_sigma"])
	numnotfixed_covsigma  = sum(poped_db["parameters"]["notfixed_covsigma"])
	
	n_fixed_eff = numnotfixed_bpop
	n_rand_eff = numnotfixed_d + numnotfixed_covd + numnotfixed_docc + numnotfixed_covdocc + numnotfixed_sigma + numnotfixed_covsigma
	fim_size = n_fixed_eff+n_rand_eff
	return fim_size