"""
#' MATLAB feval function
#' 
#' This is just a wrapper for the \code{\link{do.call}} function to behave like the feval function in MATLAB.
#' 
#' @param file.name A function or a string that is the name of a function. 
#' @param ... Arguments for the function.  Multiple arguments separated by a comma.
#' 
#' @return Output from the defined function.
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_feval.R
#' @export
#' @keywords internal
#' 
## Function written to match MATLAB function
## Author: Caiya Zhang, Yuchen Zheng
"""


def feval (func_name, *argv):
    #func_name = gsub("\\.R$","",file.name)
    if callable(func_name):
        return eval(func_name.__name__ + str(argv))
    elif type(func_name) is str:
        return eval(func_name + str(argv))




    