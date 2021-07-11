"""
## Function translated automatically using 'matlab.to.r()'


## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.cell import cell
from project.size import size
from project.zeros import zeros
from project.feval import feval
from project.LinMatrixL import LinMatrixL
from project.LinMatrixH import LinMatrixH
from project.LinMatrixLH import LinMatrixLH
from project.LinMatrixL_occ import LinMatrixL_occ
from project.hessian_eta_complex import hessian_eta_complex
from project.trace_matrix import trace_matrix
from project.diag_matlab import diag_matlab


def v(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped_db):
  # number of samples X number of samples (per individual)
    
  bUseAutoCorrelation = len(poped_db["model"]["auto_pointer"]) != 0
    
  bUseFullSigmaCorrelation = False
    
  if poped_db["settings"]["m2_switch"][1] == 0 or poped_db["settings"]["m2_switch"][1] == 1:
    returnArgs = LinMatrixL(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped_db) 
    l = returnArgs[[0]]
    poped_db = returnArgs[[1]]
    returnArgs = LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped_db) 
    h = returnArgs[[0]]
    poped_db = returnArgs[[1]]
    
    ret = zeros(0,1)
    
    if len(sigma) != 0 and bUseFullSigmaCorrelation: #Update sigma to be fully correlated
      for i in range(0, size(sigma)[0]):
        for j in range(0, size(sigma)[0]):
          if i != j:
            sigma[i,j] = np.sqrt(sigma[i,i]*sigma[j,j])
              
    #Add all IIV
    if len(d) != 0:
      ret = np.matmul(np.matmul(l, d), np.transpose(l))
    else:
      ret = zeros(len(xt_ind), len(xt_ind))
    
    if poped_db["settings"]["bUseSecondOrder"]:
      var_eta = zeros(1, len(xt_ind))
      for o in range(0,len(xt_ind)):
        hessian_eta = hessian_eta_complex(model_switch,xt_ind[o],x,a,bpop,b_ind,bocc_ind,poped_db)
        var_eta[o] = 1/4 * trace_matrix(hessian_eta * d * (2 * hessian_eta) * d)
      ret = ret + diag_matlab(var_eta)
    
    locc = cell(1, poped_db["parameters"]["NumOcc"])
    
    #Add all occasion variability
    for i in range(0, poped_db["parameters"]["NumOcc"]):
      if poped_db["parameters"]["NumOcc"] == 0:
        continue
      returnArgs = LinMatrixL_occ(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,i,poped_db) 
      locc_tmp = returnArgs[[0]]
      poped_db = returnArgs[[1]]
      if len(ret) == 0:
        ret = np.matmul(np.matmul(locc_tmp,docc), np.transpose(locc_tmp))
      else: 
        ret = ret + np.matmul(np.matmul(locc_tmp,docc), np.transpose(locc_tmp))
      
      locc[[i]] = locc_tmp
    
    if len(sigma) != 0: #If we have any residual variance
      interact = zeros(0,1) 
      if len(d) != 0: #Calculate the interaction terms
          returnArgs = LinMatrixLH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,size(sigma,1),poped_db) 
          lh = returnArgs[[0]]
          poped_db = returnArgs[[1]]
    ###numpy broadcasting???
          interact = lh
          d_sigma_prod = np.empty([d.shape[0], d.shape[1]])
          for j in range(0, d_sigma_prod.shape[1]): #frame col
            for i in range(0, d_sigma_prod.shape[0]): #frame row
              d_sigma_prod[i,j] = np.array(d[i,j]*sigma)
              i = i + 1
            j = j + 1
            i = 0
          np.matmul(np.matmul(interact, d_sigma_prod), np.transpose(lh))
          if sum(interact.size) == 2:
            ret = ret + interact
          else:
            ret = ret + diag_matlab(diag_matlab(interact))

      if bUseAutoCorrelation == True:#Add autocorrelation
        autocorr = feval(poped_db["model"]["auto_pointer"],h,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,l,locc,interact,poped_db)
        if len(ret) == 0:
          ret = autocorr
        else:
          ret = ret + autocorr
      else:#Add linearized residual model
        full_sig = np.matmul(np.matmul(h,sigma), np.transpose(h))
        if sum(full_sig.size) == 2:
          sig_tmp = full_sig
        else:
          sig_tmp = diag_matlab(diag_matlab(full_sig))
        if len(ret) == 0:
          ret = sig_tmp
        else:
          ret = ret + sig_tmp
  else:
    if poped_db["settings"]["m2_switch"][1] == 20:
      raise Exception("Analytic variance not defined")
      #ret = analytic_variance(xt_ind,x,a,bpop,d)
    else:
      raise Exception("Unknown derivative option for variance")
    
  return [ret, poped_db] 



## %This function fixes a bug in FreeMat 4.0
## function ret=diag(a)
##     if (~isempty(a) && size(a,1)==1 && size(a,2)==1)
##         ret=builtin('diag',[a]);
##     else
##         ret=builtin('diag',a);
##     end
## end

