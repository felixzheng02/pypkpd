"""


Author: Caiya Zhang, Yuchen Zheng
"""


from project.grad_all import grad_all




def gradfg_occ(x,a,bpop,b_ind,bocc_ind,currentOcc,poped_db):
#
#
# size: (number of g's x NumOccVariables)
#
# deriv of g's w$r.t. bocc's eval at b_ind, bocc_ind at occ currentOcc
#
    dfg_db0 = grad_all(poped_db["model"]["fg_pointer"],5,poped_db["parameters"]["ng"],x,a,bpop,b_ind,bocc_ind,poped_db,currentOcc = currentOcc,noPopED=True)
    return dfg_db0





