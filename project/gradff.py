"""
## Function translated automatically using 'matlab.to.r()'


## Author: Caiya Zhang, Yuchen Zheng
"""


from numpy import ndarray
from project.size import size
from project.feval import feval
from project.zeros import zeros
from project.grad_all import grad_all

def gradff(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped_db):
    #----------Model linearization with respect to random var.
    #
    # size of return is (samples per individual x number of g's)
    #
    # derivative of model w.r.t. g eval at b=b_ind
    #
    #
    fg0 = feval(poped_db["model"]["fg_pointer"],x,a,bpop,b_ind,bocc_ind)
    epsi0 = zeros(1,poped_db["parameters"]["notfixed_sigma"].size)

    dff_dg0 = grad_all(poped_db["model"]["ferror_pointer"],3,size(xt_ind,1),model_switch,xt_ind,fg0,epsi0,poped_db)

    return [dff_dg0, poped_db]

