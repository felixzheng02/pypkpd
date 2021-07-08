"""
## This python script defines main structural models needed.


#' Author: Caiya Zhang, Yuchen Zheng
"""

import numpy as np
from project.feval import do_call

#' Structural model: one-compartment, oral absorption, multiple bolus dose, parameterized using KE.
#' 
#' This is a structural model function that encodes a  model that is 
#' one-compartment, oral absorption, multiple bolus dose, parameterized using KE.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#' 
#' @param model_switch a vector of values, the same size as \code{xt}, identifying which model 
#' response should be computed for the 
#' corresponding xt value.  Used for multiple response models.
#' @param xt a vector of independent variable values (often time).
#' @param parameters A named list of parameter values.
#' @param poped_db a poped database.  This can be used to extract information that may be needed in the model file.
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.md.KE.R
#' 
#' @export


def ff_PK_1_comp_oral_md_KE(model_switch,xt,parameters,poped_db):
    ##-- Model: One comp first order absorption
    ## -- Analytic solution for both mutiple and single dosing
    list(parameters)
    y = xt
    N = np.floor(xt/parameters["TAU"]) + 1
    y = (parameters["DOSE"]*parameters["Favail"]/parameters["V"])*(parameters["KA"]/(parameters["KA"]-parameters["KE"]))*(np.exp(-parameters["KE"]*(xt-(N-1)*parameters["TAU"]))*(1-np.exp(-N*parameters["KE"]*parameters["TAU"]))/(1-np.exp(-parameters["KE"]*parameters["TAU"]))-np.exp(-parameters["KA"]*(xt-(N-1)*parameters["TAU"]))*(1-np.exp(-N*parameters["KA"]*parameters["TAU"]))/(1-np.exp(-parameters["KA"]*parameters["TAU"])))

    return [y,poped_db]




#' Structural model: one-compartment, oral absorption, multiple bolus dose, parameterized using CL.
#' 
#' This is a structural model function that encodes a  model that is 
#' one-compartment, oral absorption, multiple bolus dose, parameterized using CL.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.md.CL.R
#' 
#' @export
def ff_PK_1_comp_oral_md_CL(model_switch,xt,parameters,poped_db):
    ##-- Model: One comp first order absorption
    ## -- Analytic solution for both mutiple and single dosing
    list(parameters)
    y = xt
    N = np.floor(xt/parameters["TAU"])+1
    y = (parameters["DOSE"]*parameters["Favail"]/parameters["V"])*(parameters["KA"]/(parameters["KA"]-parameters["CL"]/parameters["V"]))*(np.exp(-parameters["CL"]/parameters["V"]*(xt-(N-1)*parameters["TAU"]))*(1-np.exp(-N*parameters["CL"]/parameters["V"]*parameters["TAU"]))/(1-np.exp(-parameters["CL"]/parameters["V"]*parameters["TAU"]))-np.exp(-parameters["KA"]*(xt-(N-1)*parameters["TAU"]))*(1-np.exp(-N*parameters["KA"]*parameters["TAU"]))/(1-np.exp(-["KA"]*parameters["TAU"])))  

    return [y,poped_db]
  

#' Structural model: one-compartment, oral absorption, single bolus dose, parameterized using KE.
#' 
#' This is a structural model function that encodes a  model that is 
#' one-compartment, oral absorption, single bolus dose, parameterized using KE.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.sd.KE.R
#' 
#' @export

## TODO: change the parameterization to be a function option
## TODO: only use md and then turn off if single dose

def ff_PK_1_comp_oral_sd_KE(model_switch,xt,parameters,poped_db):
    ##-- Model: One comp first order absorption
    list(parameters)
    y = xt
    y = (parameters["DOSE"]*parameters["Favail"]*parameters["KA"]/(parameters["V"]*(["KA"]-["KE"])))*(np.exp(-parameters["KE"]*xt)-np.exp(-parameters["KA"]*xt))

    return [y,poped_db]
  

#' Structural model: one-compartment, oral absorption, single bolus dose, parameterized using CL.
#' 
#' This is a structural model function that encodes a  model that is 
#' one-compartment, oral absorption, single bolus dose, parameterized using CL.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.sd.CL.R
#' 
#' @export
def ff_PK_1_comp_oral_sd_CL(model_switch,xt,parameters,poped_db):
    ##-- Model: One comp first order absorption
    list(parameters)
    y = xt
    y = (parameters["DOSE"]*parameters["Favail"]*parameters["KA"]/(parameters["V"]*(parameters["KA"]-parameters["CL"]/parameters["V"])))*(np.exp(-parameters["CL"]/parameters["V"]*xt)-np.exp(-parameters["KA"]*xt))
    return [y,poped_db]
    

#' Structural model: one-compartment, single bolus IV dose, parameterized using CL driving an EMAX model with a direct effect.
#' 
#' This is a structural model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PKPD.1.comp.sd.CL.emax.R
#' 
#' @export
def ff_PKPD_1_comp_sd_CL_emax(model_switch,xt,parameters,poped_db):
    list(parameters)
    y = xt
    MS = model_switch
    
    # PK model
    CONC = parameters["DOSE"]/parameters["V"]*np.exp(-parameters["CL"]/parameters["V"]*xt) 
    
    # PD model
    EFF = parameters["E0"] + CONC*parameters["EMAX"]/(parameters["EC50"] + CONC)
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return [y,poped_db]
  

#' Structural model: one-compartment, oral absorption, multiple bolus dose, 
#' parameterized using CL driving an inhibitory IMAX model with a direct effect.
#' 
#' This is a structural model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PKPD.1.comp.oral.md.CL.imax.R
#' 
#' @export
def ff_PKPD_1_comp_oral_md_CL_imax(model_switch,xt,parameters,poped_db):
    ##-- Model: One comp first order absorption + inhibitory imax
    ## -- works for both mutiple and single dosing  
    list(parameters)
      
    y = xt
    MS = model_switch
    
    # PK model
  
    returnArgs = ff_PK_1_comp_oral_md_CL(model_switch,xt,parameters,poped_db)
    CONC = returnArgs["y"]
    
    # PD model
    EFF = parameters["E0"]*(1 - CONC*["IMAX"]/(parameters["IC50"] + CONC))
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
      
    return [y,poped_db]
  

#' RUV model:  
#' Additive and Proportional.
#' 
#' This is a residual unexplained variability (RUV) model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{fError_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' @param epsi A matrix with the same number of rows as the \code{xt} vector, columns match the numbers defined in this 
#' function.
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family RUV_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.md.CL.R
#' @export
def feps_add_prop(model_switch,xt,parameters,epsi,poped_db):
    ## -- Residual Error function
    ## -- Additive + Proportional 
    returnArgs = do_call(poped_db["model"]["ff_pointer"],[model_switch,xt,parameters,poped_db]) 
    y = returnArgs[[0]]
    poped_db = returnArgs[[1]]
    y = y*(1+epsi[:,0])+epsi[:,1]
    
    return [y,poped_db]


#' RUV model:  
#' Additive .
#' 
#' This is a residual unexplained variability (RUV) model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{fError_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' @param epsi A matrix with the same number of rows as the \code{xt} vector, columns match the numbers defined in this 
#' function.
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family RUV_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_feps.add.R
#' @export
def feps_add(model_switch,xt,parameters,epsi,poped_db):
    ## -- Residual Error function
    ## -- Additive 
    returnArgs = do_call(poped_db["model"]["ff_pointer"],[model_switch,xt,parameters,poped_db]) 
    y = returnArgs[[0]]
    poped_db = returnArgs[[1]]
    y = y+epsi[:,0]
    
    return [y,poped_db]


#' RUV model:  
#' Proportional.
#' 
#' This is a residual unexplained variability (RUV) model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{fError_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' @param epsi A matrix with the same number of rows as the \code{xt} vector, columns match the numbers defined in this 
#' function.
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped_db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family RUV_models
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.sd.CL.R
#' 
#' @export
def feps_prop(model_switch,xt,parameters,epsi,poped_db):
    ## -- Residual Error function
    ## -- Proportional 
    returnArgs = do_call(poped_db["model"]["ff_pointer"],[model_switch,xt,parameters,poped_db]) 
    y = returnArgs[[0]]
    poped_db = returnArgs[[1]]
    y = y*(1+epsi[:,0])

    return [y,poped_db]


