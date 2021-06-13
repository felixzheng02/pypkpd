"""
Function translated automatically using 'matlab.to.r()'

Author: Caiya Zhang, Yuchen Zheng
"""


import project.grad_bpop as grad_bpop
import project.size as size

def m1(model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, d, sigma, poped_db):
  #
  # function computes the derivative of the
  # linerarized model function w$r.t. bpop
  # for an individual
  #
  # the output is a matrix with dimensions (ind_samps X nbpop)
  df_dbeta = grad_bpop.grad_bpop(grad_bpop.helper_LinMatrix, 5, size(xt_ind,1), model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, d, sigma, docc=None, poped_db=poped_db)

  return [df_dbeta, poped_db]
