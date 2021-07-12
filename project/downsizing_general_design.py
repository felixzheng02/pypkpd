#' Downsize a general design to a specific design
#' 
#' Function takes a design with potentially empty design 
#' variables and rescues the design so that a FIM can be calculated using \code{\link{mftot}}.
#' 
#' @param poped_db A PopED database 
#' @return A list containing:
#' \item{ni}{A vector of the number of samples in each group.}
#' \item{xt}{A matrix of sample times.  Each row is a vector of sample times for a group.}
#' \item{model_switch}{A matrix that is the same size as xt, specifying which model each sample belongs to.}
#' \item{x}{A matrix for the discrete design variables.  Each row is a group.}
#' \item{a}{A matrix of covariates.  Each row is a group.}
#' \item{bpop}{A matrix of fixed effect parameter values.}
#' 
#' @family poped_input
# @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
# @example tests/testthat/examples_fcn_doc/examples_downsizing_general_design.R
# @export
#' @keywords internal
#' 
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

import numpy as np
from project.size import size
from project.zeros import zeros

def downsizing_general_design(poped_db):
    # ------------- downsizing of general design
    
    ni = poped_db["design"]["ni"][1:poped_db["design"]["m"],:]
    xt = poped_db["design"]["xt"][1:poped_db["design"]["m"],1:max(poped_db["design_space"]["maxni"])]
    model_switch = poped_db["design"]["model_switch"][0:poped_db["design"]["m"],0:max(poped_db["design_space"]["maxni"])]
    
    if size(poped_db["design"]["x"])[1] != 0:
        x = poped_db["design"]["x"][1:poped_db["design"]["m"],0:size(poped_db["design"]["x"])[1]]
    else:
        x = zeros(poped_db["design"]["m"], 0)
    
    if size(poped_db["design"]["a"])[1] != 0:
        a = poped_db["design"]["a"][0:poped_db["design"]["m"],0:size(poped_db["design"]["a"])[1]]
        maxa = poped_db["design_space"]["maxa"][0:poped_db["design"]["m"],0:size(poped_db["design"]["a"])[1]]
        mina = poped_db["design_space"]["mina"][0:poped_db["design"]["m"],0:size(poped_db["design"]["a"])[1]]
    else:
        a = zeros(poped_db["design"]["m"],0)
        maxa = np.zeros(1)
        mina = np.zeros(1)
    
    bpop = poped_db["parameters"]["bpop"][0:poped_db["parameters"]["nbpop"],0:2]
    n = np.matmul(np.transpose(ni),np.ones(poped_db["design"]["m"],1))
    
    if poped_db["parameters"]["NumRanEff"] != 0:
        d=poped_db["parameters"]["d"][0:poped_db["parameters"]["NumRanEff"],0:2]
    else:
        d=poped_db["parameters"]["d"]
    
    maxxt=poped_db["design_space"]["maxxt"][0:poped_db["design"]["m"],0:max(poped_db["design_space"]["maxni"])]
    minxt=poped_db["design_space"]["minxt"][0:poped_db["design"]["m"],0:max(poped_db["design_space"]["maxni"])]
    
    return {"ni": ni, "xt": xt, "model_switch": model_switch, "x": x, "a": a, "bpop": bpop, 
                "n": n, "d": d, "maxxt": maxxt, "minxt": minxt, "maxa": maxa, "mina": mina}

