"""
## Result function for optimization routines
## 
## Create some output to the screen and a text file that summarizes the problem you solved.
## 
## @inheritParams RS_opt
## @inheritParams evaluate.fim
## @inheritParams Doptim
## @inheritParams create.poped.database
## @inheritParams blockexp
## @inheritParams blockheader
## @param fmf_init Initial FIM.
## @param dmf_init Initial OFV.
## @param param_cvs_init The initial design parameter RSE values in percent.
## 
## @family Helper
## @example test/test_blockfinal.R
## @export
## @keywords internal
## 
##
## Then manually adjusted to make work
## Author: Caiya Zhang, Yuchen Zheng
"""

import path
import re
import numpy as np
import warnings
from project.size import size
from matpy.matrix import Matrix
from project.tictoc import toc
from project.get_cv import get_rse
from project.ofv_criterion import ofv_criterion
from project.get_unfixed_params import get_unfixed_params

def blockfinal(fn,fmf,dmf,groupsize,ni,xt,x,a,model_switch,bpop,d,docc,sigma,poped_db,
                       fmf_init=None,dmf_init=None,param_cvs_init=None,
                       compute_inv=True,out_file=None,trflag=True,footer_flag=True,
                       run_time = None,
                       ):

    opt_xt=poped_db["settings"]["optsw"][1],opt_a=poped_db["settings"]["optsw"][3],opt_x=poped_db["settings"]["optsw"][2],
    opt_inds=poped_db["settings"]["optsw"][4],
    time_value = None
  
    if trflag is False:
        return
    if footer_flag:
        if fmf is None:
            compute_inv = False
        if type(fmf) is Matrix:
            compute_inv = False
    
    
    
    print(fn,'===============================================================================\nFINAL RESULTS\n')
    if(fn != ""):
        print('===============================================================================\nFINAL RESULTS\n')
    
    time_value = run_time

    if time_value is None:
        time_value = toc(echo=False,name=".poped_total_time")
    if opt_xt == True:
        print_xt(xt,ni,model_switch,fn,head_txt="Optimized Sampling Schedule\n")
        if fn != "":
            print_xt(xt,ni,model_switch,head_txt="\nOptimized Sampling Schedule\n")
    
    if opt_x == True:
      #     fprintf(fn,'x :\n')
      #     fprintf(fn,'%g \n',x)
      #     cat("Optimized x values:\n")
      #     print(x)
        tmp_txt = "\nOptimized Discrete Variables"
        tmp_txt = (tmp_txt + ':\n')
        print(fn, tmp_txt)

        if fn != "":
            print(tmp_txt)

        for ct1 in range(0, poped_db["design"]["m"]):
            print(fn,'Group %g: '% ct1)
            if fn != "":
                print('Group %g: '% ct1)

            for ct2 in range(0, size(poped_db["design"]["x"])):
                tmp_txt = '%g'
                if ct2 < size(poped_db["design"]["x"]):
                    tmp_txt = (tmp_txt + ' : ')        
                print(fn, tmp_txt, x[ct1,ct2])

                if fn!="":
                    print(tmp_txt, x[ct1,ct2])
            
            print(fn,'\n')
            if fn != "":
                print('\n')


    if opt_a == True:
        tmp_txt = "\nOptimized Covariates"
        tmp_txt = (tmp_txt +':\n')
        print(fn, tmp_txt)

        if fn != "":
            print(tmp_txt)

        for ct1 in range(0, poped_db["design"]["m"]):
            print(fn, 'Group %g: ', ct1)
            if fn != "":
                print('Group %g: ', ct1)
            for ct2 in range(0, size(poped_db["design"]["a"])):
                tmp_txt = '%g'
                if ct2 < size(poped_db["design"]["a"]):
                    tmp_txt = (tmp_txt + ' : ')
                print(fn, tmp_txt, a[ct1,ct2])
                if fn!="":
                    print(tmp_txt, a[ct1,ct2])
            
            print(fn,'\n')
            if fn != "":
                print('\n')
        
    
    if opt_inds == True:
        tmp_txt = "\nOptimized groupsize"
        tmp_txt = (tmp_txt + ':\n')
        print (fn, tmp_txt)
        if fn != "":
            print(tmp_txt)
        for ct1 in range(0, poped_db["design"]["m"]):
            print(fn,'Group %g: ', ct1)
            if fn != "":
                print('Group %g: ', ct1)
            for ct2 in range(0, size(poped_db["design"]["groupsize"])):
                tmp_txt = '%g'
                if ct2 < size(poped_db["design"]["groupsize"]):
                    tmp_txt = (tmp_txt + ' : ')
                print(fn, tmp_txt,groupsize[ct1,ct2])
                if fn != "":
                    print(tmp_txt,groupsize[ct1,ct2])
            print(fn,'\n')
            if fn != "":
                print('\n')
      
    
    if poped_db["settings"]["d_switch"] == True and (fn != "" or trflag > 1):
        print(fn,'\n FIM: \n')
        #write_matrix(fn,fmf)
        #MASS::write_matrix(fmf, file=fn)
        print(fn, '\n\nInverse(FIM):\n')
        #write_matrix(fn,inv(fmf))
        if compute_inv:
            inv_fmf = np.linalg.inv(fmf)
            #MASS::write_matrix(np.linalg.inv(fmf), file=fn)

    print(fn,'\nOFV = %g\n', dmf)

    if fn!="":
        print('\nOFV = %g\n',dmf)
    
    if compute_inv:
        fim=fmf_init
        param_cvs = get_rse(poped_db, bpop, np.diag(d),docc, sigma, fim=fim)
    
    output = get_unfixed_params(poped_db)
    params = output["all"]
    npar = len(params)
    
    if fn != "" or trflag > 1:
        print(fn,'\nEfficiency criterion [usually defined as det(FIM)^(1/npar)]  = %g\n',
            ofv_criterion(dmf,npar,poped_db))
    
    """
    if dmf_init is None:
        eff = efficiency(dmf_init, dmf, poped_db)
        print(fn,"\nEfficiency: \n  (%s) = %.5g\n", attr(eff,"description"),eff,both=True)
    """
    
    if param_cvs_init is None and fmf_init is not None and type(fmf_init) is Matrix and compute_inv:
        if np.isfinite(dmf_init):
            #param_vars_init=diag_matlab(inv(fmf_init))
            #returnArgs =  get_cv(param_vars_init,bpop,d,docc,sigma,poped_db) 
            #params_init = returnArgs[[1]]
            #param_cvs_init = returnArgs[[2]]
            param_cvs_init = get_rse(poped_db,bpop,np.diag(d),docc,sigma,fim=fmf_init)
        else:
            param_cvs_init = warnings.warn(get_rse(poped_db,bpop,np.diag(d),docc,sigma,fim=fmf_init))

    
    if compute_inv:
      parnam = get_parnam(poped_db)
      print(fn,'\nExpected relative standard error\n(%sRSE, rounded to nearest integer):\n','%')
      if fn != "":
        print('\nExpected relative standard error\n(%sRSE, rounded to nearest integer):\n','%')
      df = Matrix(np.array([parnam, ("%6.3g"+params), #"Variance"=param_vars, 
                       round(param_cvs_init), round(param_cvs)]),
                       (1, 4),
                       ["Parameter","Values", #"Variance"=param_vars, 
                       "RSE_0","RSE"],
                       None, None)
      #print(df,digits=3, print_gap=3,row_names=F)
      #if fn != "":
        #capture_output(print(df,digits=3, print.gap=3,row.names=F),file=fn)
    
    if time_value is not None:
      print(fn,'\nTotal running time: %g seconds\n', time_value)
      if fn != "":
        print('\nTotal running time: %g seconds\n', time_value)
    
   # end footer_flag
  #if any(class(out_file)=="file") and (fn != '')) 
    #close(fn)
  #return(invisible(time_value)) 


def print_xt(xtopt, ni, model_switch,fn="",head_txt="Optimized sample times:\n",xt_other=None,
                      digits=4):
  print(head_txt + fn)
  for j in range(0, size(xtopt)):
    xtopt_i = xtopt[j, 0:ni[j]]
    model_switch_i = model_switch[j,1:ni[j]]
    if xt_other is not None:
      xt_other_i = xt_other[j, 0:ni[j]]

    for i in np.unique(np.array([model_switch_i])):
      xtopt_i_sort = sorted(xtopt_i[model_switch_i==i])
      if xt_other is not None:
        xt_other_i_sort = xt_other_i[np.argsort(xtopt_i[model_switch_i==i])]
      # if(size(xtopt,1)>1) cat(sprintf("Group %g : ", j),file=fn)
      print("Group %g: " % j)
      if len(np.unique(np.array([model_switch_i]))) > 1:
        print("Model %g: " % i)

      if xt_other is not None:
        print((("%" + (digits + 2) + "." + digits + "g") + xt_other_i_sort) + fn)
      else:
        print((("%" + (digits + 2) + "." + digits + "g") + xtopt_i_sort) + fn)
      print("\n" + fn)
  return



def get_parnam(poped_db):
  nbpop = len(poped_db["parameters"]["notfixed_bpop"])
  nd = len(poped_db["parameters"]["notfixed_d"])
  ncovd = len(poped_db["parameters"]["notfixed_covd"])
  ndocc = len(poped_db["parameters"]["notfixed_docc"])
  ncovdocc = len(poped_db["parameters"]["notfixed_covdocc"])
  nsigma = len(poped_db["parameters"]["notfixed_sigma"])
  ncovsigma = len(poped_db["parameters"]["notfixed_covsigma"])

  not_fixed = Matrix(np.array([poped_db["parameters"]["notfixed_bpop"], 
                                poped_db["parameters"]["notfixed_d"],
                                poped_db["parameters"]["notfixed_covd"],
                                poped_db["parameters"]["notfixed_docc"],
                                poped_db["parameters"]["notfixed_covdocc"],
                                poped_db["parameters"]["notfixed_sigma"],
                                poped_db["parameters"]["notfixed_covsigma"]]),
                        (1, 7),
                        ["bpop","D","D_cov","D_occ","D_occ_cov","SIGMA","SIGMA_cov"],
                        None, None)
  
  bpop_names = poped_db["parameters"]["bpop"].get_axisnam()[0]
  d_names = poped_db["parameters"]["d"].get_axisnam()[0]
  sig_names = poped_db["parameters"]["sigma"].get_axisnam()[0]
  
  parnam = np.array([])
  for i in range(0, size(not_fixed)):
    if len(not_fixed[[i]])==0:
      next
    for j in range(0, len(not_fixed[[i]])):
      if not_fixed[i][j] == 1: 
        if (not_fixed[i]).get_axisnam() == "bpop":
          default_name = True
          if bpop_names is not None:
            if bpop_names[j] != "":
              default_name = False
              parnam = np.array([parnam, bpop_names[j]])
          if default_name:
            parnam = np.array([parnam, (not_fixed[i] +"[",j,"]")])    
        
        if any((not_fixed[i]).get_axisnam()) == "D":
          default_name = True
          if d_names is not None:
            if d_names[j] != "":
              default_name = False
              parnam = np.array([parnam, ("d_" + d_names[j])])
          if default_name:
            parnam = np.array([parnam, (not_fixed[i]).get_axisnam() +"["+j+","+j+"]"])
        
        if any((not_fixed[i]).get_axisnam()) == "SIGMA":
          default_name = True
          if sig_names is not None:
            if sig_names[j] != "":
              default_name = False
              parnam = np.array([parnam, ("sig_"+sig_names[j])])
            
          if default_name:
            parnam = np.array([parnam, (not_fixed[i]).get_axisnam()+"["+j+","+j+"]"])

        if any((not_fixed[i]).get_axisnam()) == "D_cov":
          mat_ind = np.where(poped_db["parameters"]["param_pt_val"]["d"][np.tril_indices(poped_db["parameters"]["param_pt_val"]["d"].shape[0],-1)])[j,:]
          parnam = np.array([parnam, ("D"+"["+mat_ind[0]+","+mat_ind[1]+"]")])
        
        if any((not_fixed[i]).get_axisnam()) == "D_occ":
          parnam = np.array([parnam, ((not_fixed[i]).get_axisnam()+"["+j+","+j+"]")])
            
        if len(re.search("_cov", (not_fixed[i]).get_axisnam())) != 0 and any((not_fixed[i]).get_axisnam()) =="D_cov":
          parnam = np.array([parnam, ((not_fixed[i]).get_axisnam()+"["+j+"]")]) 
  
  return parnam