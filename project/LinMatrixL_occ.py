"""
## Model linearization with respect to occasion variability parameters.
## 
## The function performs a linearization of the model with respect to the occasion  variability parameter..
## Derivative of model w.r.t. eta_occ, evaluated bocc_ind.
## 
## @inheritParams mftot
## @inheritParams LinMatrixH
## @param iCurrentOcc The current occasion.
## 
## @return A matrix of size (samples per individual x number of iovs)
##  
## @family FIM
## @export
## @keywords internal
## Function written to match MATLAB function LinMatrixL_occ()

## Author: Caiya Zhang
# """


import numpy as np
from project.gradff import gradff
from project.gradfg_occ import gradfg_occ

def LinMatrixL_occ(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,iCurrentOcc,poped_db):
    #
    # size: (samples per individual x number of iovs)
    #
    if poped_db["parameters"]["NumOcc"] == 0:
        y = 0
    else:
        returnArgs = gradff(model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, poped_db) 
        grad_ff_tmp = returnArgs[0]
        poped_db = returnArgs[1]
        y = np.matmul(grad_ff_tmp, gradfg_occ(x, a, bpop, b_ind, bocc_ind, iCurrentOcc, poped_db))
    
    return {"y": y, "poped_db": poped_db}



