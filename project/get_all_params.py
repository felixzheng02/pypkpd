"""
## Extract all model parameters from the poped database.
## 
## @param a pop_ed database
## @return A list containing:
## \item{bpop}{A vector of fixed effect parameter values.}
## \item{d}{A vector of between subject variability parameters}
## \item{covd}{A vector of the covariances of the between subject variability parameters.  Row major format of the lower triangular portion of the D (OMEGA) Matrix}
## \item{docc}{A vector of the between occasion variability (BOV) terms in the model}
## \item{covdocc}{A vector of the covariances between the BOV terms.  Row major of the lower triangular portion of the BOV Matrix. }
## \item{sigma}{A vector of the residual unexplained variances (RUV)}
## \item{covsigma}{A vector of the covariances between the RUV terms}
## \item{all}{A vector with all of the above, in the order of this list.}

## @export
## @keywords internal

## Function written to match MATLAB function get_all_params()


## Author: Caiya Zhang, Yuchen Zheng
"""

import path
import numpy as np
from project.size import size
from matpy.matrix import Matrix
from project.diag_matlab import diag_matlab
from project.util import get_dict_value
from project.data import data
from project.length import length


def get_all_params (pypkpd_db):
    #Return all params (in a vector all) with the specified order above


    bpop = Matrix(data(data(get_dict_value(pypkpd_db, "parameters", "bpop"))[:,1]))
    d = Matrix(data(data(get_dict_value(pypkpd_db, "parameters", "d"))[:,1]))
    docc = Matrix(data(data(get_dict_value(pypkpd_db, "parameters", "docc"))[:,1]))
    covd = Matrix(data(get_dict_value(pypkpd_db, "parameters", "covd")))
    sigma = Matrix(data(diag_matlab(get_dict_value(pypkpd_db, "parameters", "sigma"))))
    covsigma = Matrix(data(np.zeros([1, int(length(sigma)*(length(sigma)-1)/2)])))
    k = 1
    for i in range(0, size(get_dict_value(pypkpd_db, "parameters", "sigma"))[0]):
        for j in range(0, size(get_dict_value(pypkpd_db, "parameters", "sigma"))[1]):
            if i < j:
                covsigma[k] = data(get_dict_value(pypkpd_db, "parameters", "sigma"))[i,j]
                k = k + 1
    
    covdocc = data(get_dict_value(pypkpd_db, "parameters", "covdocc"))

    all_mat = {"bpop": bpop,
                    "d": d,
                    "covd": covd,
                    "docc": docc,
                    "covdocc": covdocc,
                    "sigma": sigma,
                    "covsigma": covsigma}

    included_in_all = []
    for value in all_mat.values():
        if length(value) != 0:
            if size(value)[1] != 1:
                included_in_all.append(np.transpose(data(value)))
            else:
                included_in_all.append(data(value))

    all = Matrix(np.concatenate(tuple(included_in_all)))

    all_mat["all"] = all

    return all_mat