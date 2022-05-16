"""
## Function written to match MATLAB function gradfg()


## Author: Caiya Zhang, Yuchen Zheng
"""


from project.grad_all import grad_all

def gradfg(x,a,bpop,b_ind,bocc_ind,poped_db):
#
#
# size: (number of g's x number of random effects)
#
# deriv of g's w.r.t. b's and bocc's eval at b_ind
#
  dfg_db0 = grad_all(poped_db["model"]["fg_pointer"],4,poped_db["parameters"]["ng"],x,a,bpop,b_ind,bocc_ind,poped_db,noPopED=True)
  return dfg_db0 


