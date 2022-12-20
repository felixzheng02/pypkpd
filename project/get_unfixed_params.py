"""
## Return all the unfixed parameters
## 
## all = vector of all unfixed params var derivative is a vector of 1 and 0, 1
## means derivative of parameter is taken w$r.t. variance otherwise w$r.t. sd If
## params is supplied then the parameter is taken from this vector instead of
## pypkpd_db
## @param pypkpd_db a PopED database.
## @param params If params is supplied then the parameters are taken from this vector. 
##   
## @return A list with the  parameters.  All unfixed parameters are also
##   returned in the "\code{all} output with the specified order 
##   (bpop,d,covd,docc,covdocc,sigma,covsigma). \code{var_derivative}  is a
##   vector of 1's or 0's, 1 means derivative of parameter is taken with respect
##   to the variance otherwise with respect to standard deviation.
## @export
## @keywords internal


Author: Caiya Zhang, Yuchen Zheng
"""



import numpy as np
from project.size import size
from matpy.matrix import Matrix
from project.diag_matlab import diag_matlab
from project.length import length
from project.util import get_dict_value
from project.data import data


def get_unfixed_params(pypkpd_db, params=None):
  
    if params is None:

        bpop = data(data(get_dict_value(pypkpd_db, "parameters", "bpop"))[:,1])
        d = data(data(get_dict_value(pypkpd_db, "parameters", "d"))[:,1])
        covd = data(get_dict_value(pypkpd_db, "parameters", "covd"))
        docc = data(data(get_dict_value(pypkpd_db, "parameters", "docc"))[:,1])
        covdocc = data(get_dict_value(pypkpd_db, "parameters", "covdocc"))
        sigma = data(diag_matlab(get_dict_value(pypkpd_db, "parameters", "sigma")))
        covsigma = data(np.zeros([1, int(length(sigma)*(length(sigma)-1)/2)]))
        k = 1
        for i in range(0, size(get_dict_value(pypkpd_db, "parameters", "sigma"))[0]):
            for j in range(0, size(get_dict_value(pypkpd_db, "parameters", "sigma"))[1]):
                if i < j:
                    covsigma[k] = data(get_dict_value(pypkpd_db, "parameters", "sigma"))[i,j]
                    k = k + 1
    else:
        nbpop = pypkpd_db["parameters"]["notfixed_bpop"].get_size()
        nd = pypkpd_db["parameters"]["notfixed_d"].get_size()
        ncovd = pypkpd_db["parameters"]["notfixed_covd"].get_size()
        ndocc = pypkpd_db["parameters"]["notfixed_docc"].get_size()
        ncovdocc = pypkpd_db["parameters"]["notfixed_covdocc"].get_size()
        nsigma = pypkpd_db["parameters"]["notfixed_sigma"].get_size()
        ncovsigma = pypkpd_db["parameters"]["notfixed_covsigma"].get_size()
        
        bpop = params[0:nbpop]
        d = params[nbpop:(nbpop+nd+1)]
        covd = params[(nbpop+nd):(nbpop+nd+ncovd+1)]
        docc = params[(nbpop+nd+ncovd):(nbpop+nd+ncovd+ndocc+1)]
        covdocc = params[(nbpop+nd+ncovd+ndocc):(nbpop+nd+ncovd+ndocc+ncovdocc+1)]
        sigma = params[(nbpop+nd+ncovd+ndocc+ncovdocc):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma+1)]
        covsigma = params[(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma+ncovsigma+1)]
    
    bpop = Matrix(bpop[data(get_dict_value(pypkpd_db, "parameters", "notfixed_bpop")) == 1])
    d = Matrix(d[data(get_dict_value(pypkpd_db, "parameters", "notfixed_d")) == 1])
    covd = Matrix(covd[data(get_dict_value(pypkpd_db, "parameters", "notfixed_covd")) == 1])
    docc = Matrix(docc[data(get_dict_value(pypkpd_db, "parameters", "notfixed_docc")) == 1])
    covdocc = Matrix(covdocc[data(get_dict_value(pypkpd_db, "parameters", "notfixed_covdocc")) == 1])
    sigma = Matrix(sigma[data(get_dict_value(pypkpd_db, "parameters", "notfixed_sigma")) == 1])
    covsigma = Matrix(covsigma[data(get_dict_value(pypkpd_db, "parameters", "notfixed_covsigma")) == 1])
    
    all = Matrix(np.concatenate((data(bpop), data(d), data(covd), data(docc), data(covdocc), data(sigma), data(covsigma)), axis=1))
    
    if pypkpd_db["settings"]["iFIMCalculationType"] != 4:
        var_derivative = Matrix(np.ones(size(all)))
    else:
        var_derivative = np.array([[np.repeat(1, bpop.size)], [np.repeat(1, d.size)], [np.repeat(1, covd.size)], [np.repeat(1, docc.size)], [np.repeat(1, covdocc.size)], [np.repeat(0, sigma.size)], [np.repeat(1,covsigma.size)]])
    
    unfixed_mat = {"bpop": bpop,
                    "d": d,
                    "covd": covd,
                    "docc": docc,
                    "covdocc": covdocc,
                    "sigma": sigma,
                    "covsigma": covsigma,
                    "all": all,
                    "var_derivative": var_derivative}
    
    return unfixed_mat

