"""
#' Compute the monte-carlo mean of a function
#'
#'Function computes the monte-carlo mean of a function by varying the parameter inputs to the function
#'
#' @param ofv_fcn A function with poped_db as the first input
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param doccdescr Matrix defining the IOV.
#' per row (row number = parameter_number) we should have:
#' \itemize{
#' \item column 1 the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
#'  3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal)
#' \item column 2  defines the mean of the variance.
#' \item column 3 defines the variance of the distribution (or length of uniform distribution).
#' }
#' @param user_distribution_pointer Function name for user defined distributions for E-family designs 
#' @return The mean of the function evaluated at different parameter values.
#' @export
#'
# @examples

#' Author: Caiya Zhang, Yuchen Zheng
"""


from project.feval import feval
from project.pargen import pargen


def mc_mean(ofv_fcn,poped_db,*argv):

    bpopdescr = poped_db["parameters"]["bpop"], 
    ddescr = poped_db["parameters"]["d"],
    doccdescr = poped_db["parameters"]["d"],
    user_distribution_pointer = poped_db["model"]["user_distribution_pointer"],
    ED_samp_size = poped_db["settings"]["ED_samp_size"],
    bLHS = poped_db["settings"]["bLHS"],
    ofv_sum = 0
    
    
    bpop_gen  =  pargen(bpopdescr,
                        user_distribution_pointer,
                        ED_samp_size,
                        bLHS,
                        None,#zeros(1,0),
                        poped_db)
    
    d_gen = pargen(ddescr,
                    user_distribution_pointer,
                    ED_samp_size,
                    bLHS,
                    None,
                    poped_db)
    
    docc_gen = pargen(doccdescr,
                        user_distribution_pointer,
                        ED_samp_size,
                        bLHS,
                        None,
                        poped_db)
    
    poped_db_tmp = poped_db

    for ct in range(0,ED_samp_size):
        poped_db_tmp["parameters"]["bpop"][:,1] = bpop_gen[ct,:]
        poped_db_tmp["parameters"]["d"][:,1] = d_gen[ct,:]
        poped_db_tmp["parameters"]["docc"][:,1] = docc_gen[ct,:]
        
        dmf_tmp = feval(ofv_fcn, poped_db_tmp, *argv)
        
        ofv_sum = ofv_sum + dmf_tmp
    
    ofv_mean = ofv_sum/poped_db["settings"]["ED_samp_size"]

    return ofv_mean
    