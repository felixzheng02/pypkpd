"""
## Evaluate the expectation of the Fisher Information Matrix (FIM) and the expectation of the OFV(FIM).
## 
## Compute the expectation of the FIM given the model, parameters, distributions of parameter uncertainty, design and methods defined in the 
## PopED database. 
## 
## @inheritParams evaluate.fim
## @inheritParams Doptim
## @inheritParams create.poped.database
## @param xtoptn The xtoptn value
## @param xoptn The xoptn
## @param aoptn The aoptn value
## @param calc_fim Should the FIM be calculated or should we just use the user defined ed_penalty_pointer.
## 
## @return A list containing the E(FIM) and E(OFV(FIM)) and the a poped_db.
## 
## @family FIM
## @family E-family
##
## @export
## @keywords internal
## 
## Function written to match MATLAB function ed_mftot()


## Author: Caiya Zhang, Yuchen Zheng
"""


from project.cell import cell
from project.zeros import zeros
from project.feval import feval
from project.pargen import pargen
from project.mftot import mftot
from project.getfulld import getfulld
from project.ofv_fim import ofv_fim

def ed_mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped_db,calc_fim=True,*argv):
    #+++++++++++++++++++++++ ED OFV(MF) VALUE
    s = 0
    s1 = 0
    
    fim_list = cell(1, poped_db["settings"]["ED_samp_size"])
    d_gen_list = cell(1, poped_db["settings"]["ED_samp_size"])
    docc_gen_list = cell(1, poped_db["settings"]["ED_samp_size"])
    
    
    bpop_gen  =  pargen(bpopdescr, poped_db["model"]["user_distribution_pointer"],
                        poped_db["settings"]["ED_samp_size"], poped_db["settings"]["bLHS"], zeros(1,0), poped_db)
    
    if calc_fim is True:
        for ct in range(0, poped_db["settings"]["ED_samp_size"]):
            d_gen = getfulld(pargen(ddescr, poped_db["model"]["user_distribution_pointer"], 1, poped_db["settings"]["bLHS"], ct, poped_db), covd)
            docc_gen = getfulld(pargen(docc, poped_db["model"]["user_distribution_pointer"], 1, poped_db["settings"]["bLHS"], ct, poped_db), poped_db["parameters"]["covdocc"])
            returnArgs =  mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop_gen[ct,],d_gen,sigma,docc_gen,poped_db,...) 
            mftmp = returnArgs[0]
            poped_db = returnArgs[1]
            s = s + ofv_fim(mftmp, poped_db)
            s1 = s1 + mftmp
            fim_list["ct"] = mftmp
            d_gen_list["ct"] = d_gen
            docc_gen_list["ct"] = docc_gen
    
    if len(poped_db["settings"]["ed_penalty_pointer"]) != 0:
        returnArgs = feval(poped_db["settings"]["ed_penalty_pointer"], fim_list, bpop_gen, d_gen_list,
                            docc_gen_list, model_switch, groupsize, ni, xtoptn, xoptn, aoptn,
                            bpopdescr, ddescr, covd, sigma, docc, poped_db) 
        ED_fim = returnArgs[0]
        ED_ofv = returnArgs[1]
        poped_db = returnArgs[2]
    else:
        ED_ofv = s/poped_db["settings"]["ED_samp_size"]
        ED_fim = s1/poped_db["settings"]["ED_samp_size"]
    
    return {"ED_fim": ED_fim, "ED_ofv": ED_ofv, "poped_db": poped_db}


