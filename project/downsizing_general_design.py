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

def downsizing_general_design(pypkpd_db):
    # ------------- downsizing of general design
    
    m = get_dict_value(pypkpd_db, "design", "m").get_value()
    ni = get_dict_value(pypkpd_db, "design", "ni").get_partial_matrix([[0, m],[None, None]])
    xt = pypkpd_db["design"]["xt"].get_partial_matrix([[0, m],[0, np.max(get_dict_value(pypkpd_db, "design_space", "maxni").get_data())]])
    xt = pypkpd_db["design"]["xt"][1:pypkpd_db["design"]["m"],1:max(pypkpd_db["design_space"]["maxni"])]
    model_switch = pypkpd_db["design"]["model_switch"][0:pypkpd_db["design"]["m"],0:max(pypkpd_db["design_space"]["maxni"])]
    
    if size(pypkpd_db["design"]["x"])[1] != 0:
        x = pypkpd_db["design"]["x"][1:pypkpd_db["design"]["m"],0:size(pypkpd_db["design"]["x"])[1]]
    else:
        x = zeros(pypkpd_db["design"]["m"], 0)
    
    if size(pypkpd_db["design"]["a"])[1] != 0:
        a = pypkpd_db["design"]["a"][0:pypkpd_db["design"]["m"],0:size(pypkpd_db["design"]["a"])[1]]
        maxa = pypkpd_db["design_space"]["maxa"][0:pypkpd_db["design"]["m"],0:size(pypkpd_db["design"]["a"])[1]]
        mina = pypkpd_db["design_space"]["mina"][0:pypkpd_db["design"]["m"],0:size(pypkpd_db["design"]["a"])[1]]
    else:
        a = zeros(pypkpd_db["design"]["m"],0)
        maxa = np.zeros(1)
        mina = np.zeros(1)
    
    bpop = pypkpd_db["parameters"]["bpop"][0:pypkpd_db["parameters"]["nbpop"],0:2]
    n = np.matmul(np.transpose(ni),np.ones(pypkpd_db["design"]["m"],1))
    
    if pypkpd_db["parameters"]["NumRanEff"] != 0:
        d=pypkpd_db["parameters"]["d"][0:pypkpd_db["parameters"]["NumRanEff"],0:2]
    else:
        d=pypkpd_db["parameters"]["d"]
    
    maxxt=pypkpd_db["design_space"]["maxxt"][0:pypkpd_db["design"]["m"],0:max(pypkpd_db["design_space"]["maxni"])]
    minxt=pypkpd_db["design_space"]["minxt"][0:pypkpd_db["design"]["m"],0:max(pypkpd_db["design_space"]["maxni"])]
    
    return {"ni": ni, "xt": xt, "model_switch": model_switch, "x": x, "a": a, "bpop": bpop, 
                "n": n, "d": d, "maxxt": maxxt, "minxt": minxt, "maxa": maxa, "mina": mina}

