'''
## Downsize a general design to a specific design
## 
## Function takes a design with potentially empty design 
## variables and rescues the design so that a FIM can be calculated using \code{\link{mftot}}.
## 
## @param pypkpd_db A PopED database 
## @return A list containing:
## \item{ni}{A vector of the number of samples in each group.}
## \item{xt}{A Matrix of sample times.  Each row is a vector of sample times for a group.}
## \item{model_switch}{A Matrix that is the same size as xt, specifying which model each sample belongs to.}
## \item{x}{A Matrix for the discrete design variables.  Each row is a group.}
## \item{a}{A Matrix of covariates.  Each row is a group.}
## \item{bpop}{A Matrix of fixed effect parameter values.}
## 
## @family poped_input
## @export
## @keywords internal
## 
## Function written to match MATLAB function downsizing_general_design()


## Author: Andrew Hooker
'''
import numpy as np
from project.size import size
from project.util import get_dict_value
from project.zeros import zeros
from matpy.matrix import Matrix

def downsizing_general_design(pypkpd_db):
    # ------------- downsizing of general design
    
    m = get_dict_value(pypkpd_db, "design", "m").get_value()
    ni = get_dict_value(pypkpd_db, "design", "ni").get_partial_matrix([[0, m],[None, None]])
    xt = get_dict_value(pypkpd_db, "design", "xt").get_partial_matrix([[0, m],[0, np.max(get_dict_value(pypkpd_db, "design_space", "maxni").get_data())]])
    model_switch = get_dict_value(pypkpd_db, "design", "model_switch").get_partial_matrix([[0, m],[0, np.max(get_dict_value(pypkpd_db, "design_space", "maxni").get_data())]])
    
    if get_dict_value(pypkpd_db, "design", "x").get_shape(1) != 0:
        x = get_dict_value(pypkpd_db, "design", "x").get_partial_matrix([[0, m], [0, None]])
    else:
        x = Matrix(np.zeros([m, 0])) # different from R: np.zeros([1, 0]) -> array([], shape=(1, 0), dtype=float64)
    
    if get_dict_value(pypkpd_db, "design", "a").get_shape(1) != 0:
        a = get_dict_value(pypkpd_db, "design", "a").get_partial_matrix([[0, m], [0, None]])
        maxa = get_dict_value(pypkpd_db, "design", "maxa").get_partial_matrix([[0, m], [0, get_dict_value(pypkpd_db, "design", "a").get_shape(1)]])
        mina = get_dict_value(pypkpd_db, "design", "mina").get_partial_matrix([[0, m], [0, get_dict_value(pypkpd_db, "design", "a").get_shape(1)]])
    else:
        a = Matrix(np.zeros(m, 0))
        maxa = Matrix(np.zeros([0, 0])) # different from R
        mina = Matrix(np.zeros([0, 0]))
    
    bpop = get_dict_value(pypkpd_db, "parameters", "bpop").get_partial_matrix([[0, get_dict_value(pypkpd_db, "parameters", "nbpop")], [0, 3]])
    n = Matrix(np.matmul(np.transpose(ni.get_data()), np.ones(m, 1)))
    
    if get_dict_value(pypkpd_db, "parameters", "NumRanEff") != 0:
        d = get_dict_value(pypkpd_db, "parameters", "d").get_partial_matrix([[0, get_dict_value(pypkpd_db, "parameters", "NumRanEff")], [0, 3]])
    else:
        d = get_dict_value(pypkpd_db, "parameters", "d")
    
    maxxt = get_dict_value(pypkpd_db, "design_space", "maxxt").get_partial_matrix([[0, m], [0, np.max(get_dict_value(pypkpd_db, "design_space", "maxni").get_data())]])
    minxt = get_dict_value(pypkpd_db, "design_space", "minxt").get_partial_matrix([[0, m], [0, np.max(get_dict_value(pypkpd_db, "design_space", "maxni").get_data())]])
    
    return {"ni": ni, "xt": xt, "model_switch": model_switch, "x": x, "a": a, "bpop": bpop, 
                "n": n, "d": d, "maxxt": maxxt, "minxt": minxt, "maxa": maxa, "mina": mina}

