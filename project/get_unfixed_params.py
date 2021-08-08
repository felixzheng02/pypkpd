"""
#' Return all the unfixed parameters
#' 
#' all = vector of all unfixed params var derivative is a vector of 1 and 0, 1
#' means derivative of parameter is taken w$r.t. variance otherwise w$r.t. sd If
#' params is supplied then the parameter is taken from this vector instead of
#' poped_db
#' @param poped_db a PopED database.
#' @param params If params is supplied then the parameters are taken from this vector. 
#'   
#' @return A list with the  parameters.  All unfixed parameters are also
#'   returned in the "\code{all} output with the specified order 
#'   (bpop,d,covd,docc,covdocc,sigma,covsigma). \code{var_derivative}  is a
#'   vector of 1's or 0's, 1 means derivative of parameter is taken with respect
#'   to the variance otherwise with respect to standard deviation.
#' @export
#' @keywords internal


Author: Caiya Zhang, Yuchen Zheng
"""



import numpy as np
from project.size import size
from project.zeros import zeros
from matpy.matrix import matrix
from project.diag_matlab import diag_matlab


def get_unfixed_params(poped_db,params=None):
  
    if params is not None:

        #type: ndarray
        bpop = poped_db["parameters"]["bpop"].get_data()[:,1]
        d = poped_db["parameters"]["d"].get_data()[:,1]
        covd = poped_db["parameters"]["covd"].get_data()
        docc = poped_db["parameters"]["docc"].get_data()[:,1]
        covdocc = poped_db["parameters"]["covdocc"].get_data()
        sigma = diag_matlab(poped_db["parameters"]["sigma"]).get_data()
        covsigma = zeros(1, (sigma.size)*((sigma.size)-1)/2).get_data()
        k = 1
        for i in range(0, size(poped_db["parameters"]["sigma"])[0]):
            for j in range(0, size(poped_db["parameters"]["sigma"])[1]):
                if i < j:
                    covsigma[k] = poped_db["parameters"]["sigma"].get_data()[i,j]
                    k = k + 1
    else:
        nbpop = poped_db["parameters"]["notfixed_bpop"].get_size()
        nd = poped_db["parameters"]["notfixed_d"].get_size()
        ncovd = poped_db["parameters"]["notfixed_covd"].get_size()
        ndocc = poped_db["parameters"]["notfixed_docc"].get_size()
        ncovdocc = poped_db["parameters"]["notfixed_covdocc"].get_size()
        nsigma = poped_db["parameters"]["notfixed_sigma"].get_size()
        ncovsigma = poped_db["parameters"]["notfixed_covsigma"].get_size()
        
        bpop = params[1:nbpop]
        d = params[(1+nbpop):(nbpop+nd)]
        covd = params[(1+nbpop+nd):(nbpop+nd+ncovd)]
        docc = params[(1+nbpop+nd+ncovd):(nbpop+nd+ncovd+ndocc)]
        covdocc = params[(1+nbpop+nd+ncovd+ndocc):(nbpop+nd+ncovd+ndocc+ncovdocc)]
        sigma = params[(1+nbpop+nd+ncovd+ndocc+ncovdocc):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma)]
        covsigma = params[(1+nbpop+nd+ncovd+ndocc+ncovdocc+nsigma):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma+ncovsigma)]
    
    bpop = bpop[poped_db["parameters"]["notfixed_bpop"] == 1]
    d = d[poped_db["parameters"]["notfixed_d"] == 1]
    covd = covd[poped_db["parameters"]["notfixed_covd"] == 1]
    docc = docc[poped_db["parameters"]["notfixed_docc"] == 1]
    covdocc = covdocc[poped_db["parameters"]["notfixed_covdocc"] == 1]
    sigma = sigma[poped_db["parameters"]["notfixed_sigma"] == 1]
    covsigma = covsigma[poped_db["parameters"]["notfixed_covsigma"] == 1]
    
    all = np.array([[bpop], [d], [covd], [docc], [covdocc], [sigma], [covsigma]])
    
    if poped_db["settings"]["iFIMCalculationType"] != 4:
        var_derivative = np.array([1,size(all)])
    else:
        var_derivative = np.array([[np.repeat(1, bpop.size)], [np.repeat(1, d.size)], [np.repeat(1, covd.size)], [np.repeat(1, docc.size)], [np.repeat(1, covdocc.size)], [np.repeat(0, sigma.size)], [np.repeat(1,covsigma.size)]])
    
    
    return {"bpop": bpop, "d": d, "covd": covd, "docc": docc, "covdocc": covdocc, "sigma": sigma, "covsigma": covsigma, "all": all, "var_derivative": var_derivative}
    

