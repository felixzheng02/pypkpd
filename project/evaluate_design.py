#' Evaluate a design
#' 
#' This function evaluates the design defined in a poped database.
#' 
#' @param poped_db A poped database
#' @param ... Extra parameters passed to \code{\link{calc_ofv_and_fim}} and \code{\link{get_rse}}
#' @return A list of elements evaluating the current design.
#' @export
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_evaluate_design.R
#' @family evaluate_design
!!
def evaluate_design(poped_db, *args):
    out = calc_ofv_and_fim(poped_db,...)
    if out["fim"] is None:
        out["rse"] = None
    else:
        out["rse"] = get_rse(out["fim"],poped_db,...)
    
    if out["rse"].keys() is not None:
        rownames(out["fim"]) = out["rse"].keys()
        colnames(out["fim"]) = out["rse"].keys()
    
    return out