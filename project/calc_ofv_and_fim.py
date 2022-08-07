"""
## Calculate the Fisher Information Matrix (FIM) and the OFV(FIM) for either point values or parameters or distributions.
## 
## This function computes the expectation of the FIM and OFV(FIM) for either point values of parameter estimates
## or parameter distributions given the model, parameters, 
## distributions of parameter uncertainty, design and methods defined in the 
## PopED database.
## 
## @inheritParams evaluate.fim
## @inheritParams Doptim
## @inheritParams create.poped.database
## @param ofv The current ofv.  If other than zero then this values is simply returned unchanged.
## @param fim The current FIM. If other than zero then this values is simply returned unchanged.
## @param use_laplace Should the Laplace method be used in calculating the expectation of the OFV?  
## @param laplace.fim Should an E(FIM) be calculated when computing the Laplace approximated E(OFV).  Typically
## the FIM does not need to be computed and, if desired,  this calculation
## is done using the standard MC integration technique, so can be slow. 
## @param evaluate_fim Should the FIM be calculated?
## 
## @return A list containing the FIM and OFV(FIM) or the E(FIM) and E(OFV(FIM)) according  to the function arguments.
## 
## @family FIM
## @family E-family
## @family evaluate_FIM
## @export


## Author: Caiya Zhang, Yuchen Zheng
"""



import re
import inspect
import numpy as np
from matpy.matrix import Matrix
from project.feval import feval
from project.ofv_fim import ofv_fim
from project.mc_mean import mc_mean
from project.getfulld import getfulld
from project.evaluate_fim import evaluate_fim
from project.evaluate_e_ofv_fim import evaluate_e_ofv_fim
from project.param_set import param_set
from project.util import default_if_none
from project.util import get_dict_value


def calc_ofv_and_fim(pypkpd_db,
                        ofv = 0,
                        fim = 0, 
                        d_switch = None,  
                        bpopdescr = None, 
                        ddescr = None,
                        bpop = None, 
                        d = None, 
                        docc_full = None, 
                        model_switch = None, 
                        ni = None, 
                        xt = None, 
                        x = None, 
                        a = None, 
                        fim_calc_type = None, 
                        use_laplace = None, 
                        laplace_fim = False, 
                        ofv_fun = None,
                        evaluate_FIM = True,
                        **kwargs):
  
    d_switch = param_set(d_switch, pypkpd_db, "settings", "d_switch")  
    bpopdescr = param_set(bpopdescr, pypkpd_db, "parameters", "bpop") 
    ddescr = param_set(ddescr, pypkpd_db, "parameters", "d")
    bpop = default_if_none(bpop, bpopdescr.get_partial_matrix([[None, None], [1, 1]]))
    d = default_if_none(d, getfulld(ddescr.get_data()[:,1],pypkpd_db["parameters"]["covd"]))
    docc_full = default_if_none(docc_full, getfulld(pypkpd_db["parameters"]["docc"].get_data()[:,1], pypkpd_db["parameters"]["covdocc"]))
    model_switch = param_set(model_switch, pypkpd_db, "design", "model_switch")
    ni = param_set(ni, pypkpd_db, "design", "ni")
    xt = param_set(xt, pypkpd_db, "design", "xt")
    x = param_set(x, pypkpd_db, "design", "x")
    a = param_set(a, pypkpd_db, "design", "a")
    fim_calc_type = param_set(fim_calc_type, pypkpd_db, "settings", "iFIMCalculationType")
    use_laplace = param_set(use_laplace, pypkpd_db, "settings", "iEDCalculationType")
    ofv_fun = param_set(ofv_fun, pypkpd_db, "settings", "ofv_fun")

    ## compute the OFV
    if ofv == 0:
        if d_switch: 
            if ofv_fun is None:
                if type(fim) is not Matrix: 
                    fmf = evaluate_fim(pypkpd_db,
                                        bpop_val=bpop,
                                        d_full=d,
                                        docc_full=docc_full,
                                        sigma_full=get_dict_value(pypkpd_db, "parameters", "sigma"),
                                        model_switch=model_switch,
                                        ni=ni,
                                        xt=xt,
                                        x=x,
                                        a=a,
                                        groupsize=get_dict_value(pypkpd_db, "design", "groupsize"),
                                        fim_calc_type=fim_calc_type,
                                        **kwargs)
                
                #     returnArgs =  mftot(model_switch,pypkpd_db["design"]groupsize,ni,xt,x,a,bpop,d,pypkpd_db["parameters"]sigma,docc_full,pypkpd_db) 
                #     fmf = returnArgs[1]
                #     pypkpd_db = returnArgs[[2]]
                
                dmf = ofv_fim(fmf, pypkpd_db, **kwargs)
            else:
                ## update pypkpd_db with options supplied in function
                called_args = match_call()
                default_args = formals()
                for i in called_args.keys()[-1]:
                    if len(re.match("^poped\\.db\\$", inspect.getsource(default_args[i]))) == 1:
                        #eval(parse(text=paste(capture.output(default_args[[i]]),"=",called_args[[i]])))
                        if eval(i) is not None:
                            eval(inspect.getsource(default_args[i]) + "=" + i)
                out_tmp = eval(str(ofv_fun) + str(pypkpd_db) + str(*argv))
                dmf = out_tmp[0]
                fmf = None
                if len(out_tmp) > 1:
                    fmf = out_tmp[1]
            
        else:   # e-family
            if ofv_fun is None:
                output = evaluate_e_ofv_fim(pypkpd_db,
                                        fim_calc_type=fim_calc_type,
                                        bpop=bpopdescr,
                                        d=ddescr,
                                        covd=pypkpd_db["parameters"]["covd"],
                                        docc=pypkpd_db["parameters"]["docc"],
                                        sigma=pypkpd_db["parameters"]["sigma"],
                                        model_switch=model_switch,
                                        ni=ni,
                                        xt=xt,
                                        x=x,
                                        a=a,
                                        groupsize=pypkpd_db["design"]["groupsize"],
                                        use_laplace=use_laplace,
                                        laplace_fim=laplace_fim, 
                                        *argv)
                dmf = output["E_ofv"]
                fmf = output["E_fim"] 
            else:
                ## update pypkpd_db with options supplied in function
                called_args = match_call()
                default_args = formals()
                for i in called_args.keys()[-1]:
                    if len(re.match("^poped\\.db\\$", inspect.getsource(default_args[i]))) == 1:
                        #eval(parse(text=paste(capture.output(default_args[[i]]),"=",called_args[[i]])))
                        if eval(i) is not None: 
                            eval(inspect.getsource(default_args[i]) + "=" + i)
            
            
                dmf = mc_mean(ofv_fun, pypkpd_db, *argv)
                fmf = None
    
        ofv = dmf
        if type(fim) is Matrix:
            fim = fmf
    
    
    ## Should we compute the FIM?
    calc_fim = False

    if fim is None:
        calc_fim = True
    elif type(fim) is np.ndarray: 
        calc_fim = True
     
    if se_laplace is True and laplace_fim is False:
        calc_fim = False
    if evaluate_fim is False:
        calc_fim = False
    
    if calc_fim is True:
        if d_switch is True:
            fmf = evaluate_fim(pypkpd_db, *argv)
            bpop_val = bpop
            d_full = d
            docc_full = docc_full
            sigma_full = pypkpd_db["parameters"]["sigma"]
            model_switch = model_switch
            ni = ni
            xt = xt
            x = x
            a = a
            groupsize = pypkpd_db["design"]["groupsize"]
            fim_calc_type = fim_calc_type
        #     returnArgs =  mftot(model_switch,pypkpd_db["design"]groupsize,ni,xt,x,a,bpop,d,pypkpd_db["parameters"]sigma,docc_full,pypkpd_db) 
        #     fmf = returnArgs[1]
        #     pypkpd_db = returnArgs[[2]]
        else:
            output = evaluate_e_ofv_fim(pypkpd_db, *argv)  
            fim_calc_type = fim_calc_type
            bpop = bpopdescr
            d = ddescr
            covd = pypkpd_db["parameters"]["covd"]
            docc = pypkpd_db["parameters"]["docc"]
            sigma = pypkpd_db["parameters"]["sigma"]
            model_switch = model_switch
            ni = ni
            xt = xt
            x = x
            a = a
            groupsize = pypkpd_db["design"]["groupsize"]
            use_laplace = False
            fmf = output["E_fim"]
        fim = fmf
    
    calc_mat = Matrix(np.array([[ofv], [fim]]))
    calc_mat.set_datanam(["ofv", "fim"])
    return calc_mat
