"""
## Build poped parameter function from a model function
##
## @param model A string of text describing the model function name
## @param covariates A list of covariate names to be filtered out of the model 
## @param par_names A list of parameter names in the model file.  If not supplied then 
## all undefined variables in the model file are extracted and the covariate names are
## filtered out of that list.
## @param etas Can be "exp", "prop", "add" or "none".  Either one value for all parameters or
## a list defining the model per parameter.
## @param no_etas Parameters that should not have etas associated with them.
## @param env The environment to create the function in.
##
## @return A parameter model function to be used as input to PopED calculations.
## @export
## @importFrom codetools findGlobals
##
## @examples
## build_sfg(model="ff.PK.1.comp.oral.md.CL")
## 
## etas = c(Favail="exp",KA="exp",V="add",CL="exp")
## build_sfg(model="ff.PK.1.comp.oral.md.CL",etas = etas)


## Author: Caiya Zhang, Yuchen Zheng
"""


import re
import builtins
import numpy as np
from project.models import ff_PK_1_comp_oral_sd_CL

## create an empty function
def sfg_tmp(x1,a,bpop,b,bocc):
    return None

def build_sfg(model=ff_PK_1_comp_oral_sd_CL,
            covariates=["dose","tau"],
            par_names=None,
            etas="exp", # can be exp, prop, add, none. can be one for all or per parameter
            no_etas=["F","Favail"]):
            #env = parent_frame()):

#> <bytecode: 0x7fe20a979808>
#> <environment: namespace:PopED>
## -- parameter definition function 
## -- names match parameters in function ff

    ## get variable of function
    parameter_names_ff = par_names
    if parameter_names_ff is None: 
        parameter_names_ff = locals(model)["variables"]  
    
    ## create an empty function
    sfg_tmp()

    ## get body of function
    parameters = "parameters = np.array(["
    bpop_num = 1
    b_num = 1
    a_num = 1
    
    cov_locs = re.search(("^" + "|".joint(covariates) +'["'), parameter_names_ff, re.I)
    covariate_names = parameter_names_ff[cov_locs]
    parameter_names = parameter_names_ff
    if len(cov_locs) > 0:
        parameter_names = parameter_names_ff[-cov_locs]
    
    # match names
    df = None
    if etas is not None:
        if all(parameter_names in etas[0].keys()):
            eta_mod = etas[parameter_names]
            df = np.concatenate([parameter_names],[eta_mod],axis=1)
        
    
    if df is not None:
        eta_mod = etas
        df = np.concatenate(parameter_names,eta_mod,axis=1)
    
    if no_etas is not None or len(etas)==1:
        df[parameter_names in no_etas,"eta_mod"] = "None"
    
    
    for k in range(0,df.shape[0]):
        ending = ", "
        if k == df.shape[0] and len(covariate_names) == 0:
            ending = "])"
        
        parameters = parameters + df[k,"parameter_names"] + "=bpop[" + bpop_num + "]"
        bpop_num = bpop_num + 1
        if df[k,"eta_mod"] == "exp": 
            parameters = parameters + '*exp(b["' + b_num + '"])'
        if df[k,"eta_mod"] == "prop": 
            parameters = parameters + '*(1+b["' + b_num + '"])'
        if df[k,"eta_mod"] == "add":
            parameters = parameters + '+ b["' + b_num + '"]'  
        if df[k,"eta_mod"] == "None":
            parameters = parameters  
        if df[k,"eta_mod"] != "None":
            b_num = b_num + 1
        parameters = parameters + ending        
    
    
    if (len(covariate_names) != 0):
        for k in range(0, len(covariate_names)):
            ending = ", "
            if k == len(covariate_names): 
                ending = ")"
            parameters = parameters + covariate_names[k] + '=a["' + a_num + '"]' + ending
            a_num = a_num + 1
        
    
        
    ## add body of funciton and set environment
    '''
    text = parameters
    body(sfg_tmp) = text
    environment(sfg_tmp) = env
    '''

    return sfg_tmp