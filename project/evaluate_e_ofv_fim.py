"""
## Evaluate the expectation of the Fisher Information Matrix (FIM) and the expectation of the OFV(FIM).
## 
## Compute the expectation of the FIM and OFV(FIM) given the model, parameters, distributions of parameter uncertainty, design and methods defined in the 
## PopED database. Some of the arguments coming from the PopED database can be overwritten;  
## by default these arguments are \code{None} in the 
## function, if they are supplied then they are used instead of the arguments from the PopED database.
## 
## @inheritParams evaluate.fim
## @inheritParams Doptim
## @inheritParams create.poped.database
## @param use_laplace Should the Laplace method be used in calculating the expectation of the OFV?  
## @param laplace.fim Should an E(FIM) be calculated when computing the Laplace approximated E(OFV).  Typically
## the FIM does not need to be computed and, if desired,  this calculation
## is done using the standard MC integration technique, so can be slow. 
## 
## @return A list containing the E(FIM) and E(OFV(FIM)) and the a poped_db updated according  to the function arguments.
## 
## @family FIM
## @family E-family
## @family evaluate_FIM
## 
## @export


## Author: Caiya Zhang, Yuchen Zheng
"""


import re
import inspect
from project.ed_mftot import ed_mftot
from project.ed_laplace_ofv import ed_laplace_ofv
from project.downsizing_general_design import downsizing_general_design

def evaluate_e_ofv_fim(poped_db,*argv):
    fim_calc_type = None
    bpop = poped_db["parameters"]["bpop"]
    d = poped_db["parameters"]["d"]
    covd = poped_db["parameters"]["covd"]
    docc = poped_db["parameters"]["docc"]
    sigma = poped_db["parameters"]["sigma"]
    model_switch = None
    ni = None
    xt = None
    x = None
    a = None
    groupsize = poped_db["design"]["groupsize"]
    deriv_type = None
    bLHS = poped_db["settings"]["bLHS"],
    ofv_calc_type = poped_db["settings"]["ofv_calc_type"]
    ED_samp_size = poped_db["settings"]["ED_samp_size"]
    use_laplace = poped_db["settings"]["iEDCalculationType"]
    laplace_fim = False

    ## update poped_db with options supplied in function
    called_args = match_call()
    default_args = formals()
    for i in called_args.keys()[-1]:
        if len(re.findall('^poped\\_db\\["', inspect.getsource(default_args[[i]]))) == 1:
            #eval(parse(text=paste(capture.output(default_args[[i]]),"=",called_args[[i]])))
            eval(str(inspect.getsource(default_args[[i]])) + "=" + str(i))
        
       
    
    downsize_list = downsizing_general_design(poped_db)
    if ni is None: 
        ni = downsize_list["ni"]
    if xt is None: 
        xt = downsize_list["xt"]
    if model_switch is None:
        model_switch = downsize_list["model_switch"]
    if x is None:
        x = downsize_list["x"]
    if a is None:
        a = downsize_list["a"]    
    
    if groupsize is None:
        groupsize = poped_db["design"]["groupsize"]
    
    if fim_calc_type is not None: 
        poped_db["settings"]["iFIMCalculationType"] = fim_calc_type
    
    if deriv_type is not None:
        poped_db["settings"]["m1_switch"] = deriv_type
        poped_db["settings"]["m2_switch"] = deriv_type
        poped_db["settings"]["hle_switch"] = deriv_type
        poped_db["settings"]["gradff_switch"] = deriv_type
        poped_db["settings"]["gradfg_switch"] = deriv_type
    
    
    E_fim = None
    E_ofv = None
    
    if use_laplace is False:
        output = ed_mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,covd,sigma,docc,poped_db,*argv)
        E_fim = output["ED_fim"]
        E_ofv = output["ED_ofv"]
        poped_db=output["poped_db"]
    else:
        E_ofv  = ed_laplace_ofv(model_switch,groupsize,ni,xt,x,a,bpop,d,covd,sigma,docc,poped_db,*argv)["f"] 
        if laplace_fim is True:
            E_fim = ed_mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,covd,sigma,docc,poped_db,*argv)["ED_fim"]

    return {"E_ofv": E_ofv, "E_fim": E_fim, "poped_db": poped_db}
