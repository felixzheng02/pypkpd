"""
#' Calculate the Fisher Information Matrix (FIM) and the OFV(FIM) for either point values or parameters or distributions.
#' 
#' This function computes the expectation of the FIM and OFV(FIM) for either point values of parameter estimates
#' or parameter distributions given the model, parameters, 
#' distributions of parameter uncertainty, design and methods defined in the 
#' PopED database.
#' 
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param ofv The current ofv.  If other than zero then this values is simply returned unchanged.
#' @param fim The current FIM. If other than zero then this values is simply returned unchanged.
#' @param use_laplace Should the Laplace method be used in calculating the expectation of the OFV?  
#' @param laplace.fim Should an E(FIM) be calculated when computing the Laplace approximated E(OFV).  Typically
#' the FIM does not need to be computed and, if desired,  this calculation
#' is done using the standard MC integration technique, so can be slow. 
#' @param evaluate_fim Should the FIM be calculated?
#' 
#' @return A list containing the FIM and OFV(FIM) or the E(FIM) and E(OFV(FIM)) according  to the function arguments.
#' 
#' @family FIM
#' @family E-family
#' @family evaluate_FIM
#'  
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_ed.R
#' @example tests/testthat/examples_fcn_doc/examples_calc_ofv_and_fim.R
#' @export

#' Author: Caiya Zhang, Yuchen Zheng
"""



import re
import numpy as np
from project.feval import do_call
from project.ofv_fim import ofv_fim
from project.mc_mean import mc_mean
from project.getfulld import getfulld
from project.evaluate_e_ofv_fim import evaluate_e_ofv_fim


def calc_ofv_and_fim(poped_db, *args):
  
    ofv = 0,
    fim = 0, 
    d_switch = poped_db["settings"]["d_switch"],  
    bpopdescr = poped_db["parameters"]["bpop"], 
    ddescr = poped_db["parameters"]["d"],
    bpop = bpopdescr[:,1], 
    d = getfulld(ddescr[:,1],poped_db["parameters"]["covd"]), 
    docc_full = getfulld(poped_db["parameters"]["docc"][:,1],poped_db["parameters"]["covdocc"]), 
    model_switch = poped_db["design"]["model_switch"], 
    ni = poped_db["design"]["ni"], 
    xt = poped_db["design"]["xt"], 
    x = poped_db["design"]["x"], 
    a = poped_db["design"]["a"], 
    fim_calc_type = poped_db["settings"]["iFIMCalculationType"], 
    use_laplace = poped_db["settings"]["iEDCalculationType"], 
    laplace_fim = False, 
    ofv_fun = poped_db["settings"]["ofv_fun"],
    evaluate_fim = True,

    ## compute the OFV
    if ofv == 0:
        if d_switch: 
            if ofv_fun is not None:
                if type(fim) is np.ndarray: 
                    fmf = evaluate_fim(poped_db,*args)
                    bpop_val = bpop
                    d_full = d
                    docc_full = docc_full
                    sigma_full = poped_db["parameters"]["sigma"]
                    model_switch = model_switch
                    ni=ni
                    xt=xt
                    x=x
                    a=a
                    groupsize=poped_db["design"]["groupsize"]
                    fim_calc_type=fim_calc_type
                
                #     returnArgs =  mftot(model_switch,poped_db["design"]groupsize,ni,xt,x,a,bpop,d,poped_db["parameters"]sigma,docc_full,poped_db) 
                #     fmf = returnArgs[[1]]
                #     poped_db = returnArgs[[2]]
                
                dmf = ofv_fim(fmf,poped_db,*args)
            else:
                ## update poped_db with options supplied in function
                called_args = match_call()
                default_args = formals()
                for i in called_args.keys()[-1]:
                    if len(re.match("^poped\\.db\\$",capture.output(default_args[[i]])))==1:
                        #eval(parse(text=paste(capture.output(default_args[[i]]),"=",called_args[[i]])))
                        if eval(parse(text=paste(i))) is not None:
                            eval(parse(text=paste(capture.output(default_args[[i]]),"=",i)))
                out_tmp = do_call(ofv_fun,[poped_db,*args])
                dmf = out_tmp[[0]]
                fmf = None
                if len(out_tmp) > 1:
                    fmf = out_tmp[[1]]
            
        else:   # e-family
            if ofv_fun is None:
                output = evaluate_e_ofv_fim(poped_db,
                                        fim_calc_type=fim_calc_type,
                                        bpop=bpopdescr,
                                        d=ddescr,
                                        covd=poped_db["parameters"]["covd"],
                                        docc=poped_db["parameters"]["docc"],
                                        sigma=poped_db["parameters"]["sigma"],
                                        model_switch=model_switch,
                                        ni=ni,
                                        xt=xt,
                                        x=x,
                                        a=a,
                                        groupsize=poped_db["design"]["groupsize"],
                                        use_laplace=use_laplace,
                                        laplace_fim=laplace_fim, 
                                        *args)
                dmf = output["E_ofv"]
                fmf = output["E_fim"] 
            else:
                ## update poped_db with options supplied in function
                called_args = match.call()
                default_args = formals()
                for i in names(called_args)[-1]:
                    if(len(re.match("^poped\\.db\\$",capture_output(default_args[[i]])))==1):
                        #eval(parse(text=paste(capture.output(default_args[[i]]),"=",called_args[[i]])))
                        if eval(parse(text=paste(i))) is not None: 
                            eval(parse(text=paste(capture_output(default_args[[i]]),"=",i)))
            
            
                dmf = mc_mean(ofv_fun,poped_db,...)
                fmf = None
    
        ofv = dmf
        if type(fim) is np.np.ndarray:
            fim = fmf
    
    
    ## Should we compute the FIM?
    calc_fim = False

    if fim is None:
        calc_fim = True
    elif type(fim) is np.ndarray: 
        calc_fim = True
     
    if se_laplace is True and laplace.fim is False:
        calc_fim = False
    if evaluate_fim is False:
        calc_fim = False
    
    if calc_fim is True:
        if d_switch is True:
            fmf = evaluate_fim(poped_db,*args)
            bpop_val = bpop
            d_full = d
            docc_full = docc_full
            sigma_full = poped_db["parameters"]["sigma"]
            model_switch = model_switch
            ni = ni
            xt = xt
            x = x
            a = a
            groupsize = poped_db["design"]["groupsize"]
            fim_calc_type = fim_calc_type
        #     returnArgs =  mftot(model_switch,poped_db["design"]groupsize,ni,xt,x,a,bpop,d,poped_db["parameters"]sigma,docc_full,poped_db) 
        #     fmf = returnArgs[[1]]
        #     poped_db = returnArgs[[2]]
        else:
            output = evaluate_e_ofv_fim(poped_db,*args)  
            fim_calc_type = fim_calc_type
            bpop = bpopdescr
            d = ddescr
            covd = poped_db["parameters"]["covd"]
            docc = poped_db["parameters"]["docc"]
            sigma = poped_db["parameters"]["sigma"]
            model_switch = model_switch
            ni = ni
            xt = xt
            x = x
            a = a
            groupsize = poped_db["design"]["groupsize"]
            use_laplace = False
            fmf = output["E_fim"]
        fim = fmf
    
    return {"ofv": ofv, "fim": fim}
