"""
#' Compute efficiency.
#' 
#' Efficiency calculation between two designs.
#' 
#' 
#' @param ofv_init An initial objective function
#' @param ofv_final A final objective function.
#' @param npar The number of parameters to use for normalization.
#' @param poped_db a poped database
#' @param use_log Are the `ofv` arguments in the log space?
#' @inheritParams ofv_fim
#' @inheritParams poped_optim
#' @inheritParams create.poped.database
#' 
#' @return The specified efficiency value depending on the ofv_calc_type.  
#' The attribute "description" tells you how the calculation was made 
#' \code{attr(return_vale,"description")}
#' 
#' @family FIM
#' 
#' 
## @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
## @example tests/testthat/examples_fcn_doc/examples_ofv_criterion.R
#' 
#' @export
# Author: Caiya Zhang, Yuchen Zheng
# """


import math
from project.get_fim_size import get_fim_size


def efficiency(ofv_init,
               ofv_final,
               poped_db,
               use_log=True):
    npar = get_fim_size(poped_db)
    ofv_calc_type = poped_db["settings"]["ofv_calc_type"]
    ds_index = poped_db["parameters"]["ds_index"]
    
    eff = ofv_final/ofv_init
    # attr(eff,"description") <- "ofv_final / ofv_init"
    if ofv_calc_type == 1: # D-Optimal Design
        eff = eff^(1/npar)
        # attr(eff,"description") <- "(ofv_final / ofv_init)^(1/n_parameters)"
    if ofv_calc_type == 4: # lnD-Optimal Design
        eff = math.exp(ofv_final)/math.exp(ofv_init)^(1/npar)
        # attr(eff,"description") <- "(exp(ofv_final) / exp(ofv_init))^(1/n_parameters)"
    if ofv_calc_type == 6: # Ds-Optimal design
        if use_log:
            eff = math.exp(ofv_final)/math.exp(ofv_init)^(1/)
            # attr(eff,"description") <- "(exp(ofv_final) / exp(ofv_init))^(1/sum(interesting_parameters))"
        else:
            eff = eff^(1/)
            # attr(eff,"description") <- "(ofv_final / ofv_init)^(1/sum(interesting_parameters))"
    return eff   