"""

Author: Caiya Zhang, Yuchen Zheng
"""

import itertools
import numpy as np
import pandas as pd
from numpy.core.records import array
from scipy import stats
from project.v import v
from project.feval import feval
from project.zeros import zeros
from project.mf_all import mf_all
from project.get_fim_size import get_fim_size




def mf_all_loq(model_switch_i,xt_i,x_i,a_i,bpop_val,d_full,sigma_full,docc_full,
                       poped_db,
                       loq = -np.Infinity, # vector of length number of models
                       loq_method = 1,#poped_db["settings"]loq_method,
                       loq_PI_conf_level = 0.95,#poped_db["settings"]loq_PI_conf_level,
                       loq_prob_limit = 0.001,#poped_db["settings"]loq_prob_limit,
                       loq_start_time = None,
                       uloq = np.Infinity,
                       uloq_method = 1,
                       uloq_start_time = None,
                       verbose = False,
                       *argv):

  
    # TODO: add to poped_db
    
    # PRED calculations based on FO
    b_ind = poped_db["parameters"]["b_global"][:,1]*0
    bocc_ind = poped_db["parameters"]["bocc_global"][0]*0
    g0 = feval(poped_db["model"]["fg_pointer"], x_i, a_i, bpop_val, b_ind, bocc_ind)
    pred = feval(poped_db["model"]["ff_pointer"], model_switch_i, xt_i, g0, poped_db)
    pred = pred[0]
    pred = pred[1].np.flatten()
    
    fim_size = get_fim_size(poped_db)

    
    n_mod_tmp, idx = np.unique(np.array([model_switch_i]), return_index=True)
    n_mod = n_mod_tmp[idx.argsort()]
    
    loq_full = np.repeat(np.nan, pred.size)
    uloq_full = np.repeat(np.nan, pred.size)
    
    if loq.size == 1:
        loq_full = np.repeat(loq,pred.size)
    if uloq.size == 1: 
        uloq_full = np.repeat(uloq,pred.size)
    
    if loq.size == n_mod:
        for k in n_mod_tmp[idx.argsort()]:
            loq_full[model_switch_i == k] = loq[k]
    
    if uloq.size == n_mod:
        for k in n_mod_tmp[idx.argsort()]:
            uloq_full[model_switch_i == k] = uloq[k]
    
    if loq_start_time is not None:
        loq_full[xt_i < loq_start_time] = -np.Infinity
    if uloq_start_time is not None:
        uloq_full[xt_i < uloq_start_time] = np.Infinity
    
    if any(np.isnan(loq_full)) or any(np.isnan(uloq_full)):
        raise Exception("loq or uloq not specified properly") 
    
    # D2 method
    if uloq_method == 2:
        uloq_obs_master = pred > uloq_full
    if loq_method == 2:
        bloq_obs_master = pred < loq_full
    
    # D6 method
    if loq_method == 1 or uloq_method == 1: 
        
        # COV calculations based on FO
        cov = v(model_switch_i,xt_i,x_i,a_i,bpop_val,b_ind,bocc_ind,d_full,sigma_full,docc_full,poped_db)[0]
        
        # compute points that have PI that overlaps LOQ
        PI_alpha = 1-loq_PI_conf_level
        z_val = stats.norm.ppf(1-PI_alpha/2)
        se_val = np.sqrt(np.diag(cov))
        ci_u = pred + z_val*se_val 
        ci_l = pred - z_val*se_val 
        
        # df = tibble::tibble(pred=c(pred),ci_l=c(ci_l),ci_u=c(ci_u),loq=loq)
        # df = df %>% dplyr::mutate(above=dplyr::if_else(ci_l>loq,1,0)) %>% 
        #   dplyr::mutate(below=dplyr::if_else(ci_u<loq,1,0)) %>% 
        #   dplyr::mutate(overlap=dplyr::if_else(ci_u>loq & ci_l<loq,1,0)) %>% 
        #   dplyr::mutate(bloq_obs=dplyr::if_else(below==1 & overlap==0,1,dplyr::if_else(overlap==1,2,0)))
        
        if loq_method == 1:
            overlap = loq_full*0
            below = overlap
            above = below
        
            above[ci_l > loq_full] = 1
            below[ci_u < loq_full] = 1
            overlap[ci_u > loq_full] = 1
            overlap[ci_l < loq_full] = 1
            
            bloq_obs_master = 0*above + 2
            bloq_obs_master[below == 1 and overlap == 0] = 1
            bloq_obs_master[above == 1 and overlap == 0] = 0
        

        if uloq_method == 1:
            overlap_u = uloq_full*0
            below_u = overlap_u
            above_u = below_u

            above_u[ci_l > uloq_full] = 1
            below_u[ci_u < uloq_full] = 1
            overlap_u[ci_u > uloq_full] = 1
            overlap_u[ci_l < uloq_full] = 1
            
            uloq_obs_master = 0*above_u + 2
            uloq_obs_master[below_u == 1 and overlap_u == 0] = 0
            uloq_obs_master[above_u == 1 and overlap_u == 0] = 1
    
    #bloq_obs_master = df$bloq_obs
    #bloq_obs_master = bloq_obs_master*0+2
    
    loq_obs_master = bloq_obs_master
    loq_obs_master[uloq_obs_master == 1] = 1
    loq_obs_master[uloq_obs_master == 2 and bloq_obs_master != 1] = 2
    
    # number of potential loq_obs
    n_pot_loq = sum(loq_obs_master == 2)
    
    if n_pot_loq > 0: # D6 Method
        
        # combination of potential obs with datapoints above LOQ or below ULOQ
        s = [None]*(2*n_pot_loq)
        for i in range(0, n_pot_loq=1):
            s[i] = 0
            i = i + 1
        for j in range(n_pot_loq, 2*n_pot_loq):
            s[j] = 1
            j = j + 1
        loq_obs_init = sorted(set(itertools.permutations(s, n_pot_loq)))
        loq_obs_init = np.asarray(list(loq_obs_init))
        
        # map for type of observation
        # 0 = normal observations
        # 1 = BLOQ
        # 2 = ULOQ
        # 3 = Could be either BLOQ or ULOQ (and need expanding)
        loq_obs_map = loq_obs_master[loq_obs_master == 2]*0 + 1
        uloq_obs_map = uloq_obs_master[loq_obs_master == 2]
        bloq_obs_map = bloq_obs_master[loq_obs_master == 2]
        loq_obs_map[uloq_obs_map == 2 and bloq_obs_map != 2] = 2 
        loq_obs_map[uloq_obs_map == 2 and bloq_obs_map == 2] = 3 
        
        loq_obs_short = np.array(loq_obs_map).reshape[loq_obs_init.shape[0], loq_obs_map.size]
        loq_obs_short[loq_obs_init == 0] = 0
        
        
        # expand rows that could be BLOQ or ULOQ 
        exp_rows = np.empty(loq_obs_short.shape[0])
        for i in range(0, loq_obs_short.shape[0]):
                exp_rows[i] = any(loq_obs_short[i,:] == 3)
                i = i + 1
        loq_obs = loq_obs_short[exp_rows is False,:]
        if any(exp_rows is True):
            loq_obs_tmp = loq_obs_short[exp_rows,]
            
            # expand rows
            for i in range(0, loq_obs_tmp.shape[0]):
                #i = 1
                obs_tmp = loq_obs_tmp[i,:]

                s = [None]*(2*sum(int(obs_tmp==3)))
                for i in range(0, sum(int(obs_tmp==3))+1):
                    s[i] = 1
                    i = i + 1
                for j in range(sum(int(obs_tmp==3)), 2*sum(int(obs_tmp==3))):
                    s[j] = 2
                    j = j + 1
                perm_tmp = sorted(set(itertools.permutations(s, sum(int(obs_tmp==3)))))
                perm_tmp=np.asarray(list(perm_tmp))

                obs_tmp_exp = np.array(obs_tmp).reshape(perm_tmp.shape[0], obs_tmp.shape[1])
                obs_tmp_exp[obs_tmp_exp==3] = perm_tmp[:,:]
                loq_obs = np.concatenate(loq_obs,obs_tmp_exp,axis=0)
            
        
        # make sure that mapped values are all accounted for 
        if any(loq_obs == 3): 
            raise Exception("Combinations not fully expanded")
        
        # cat(loq_obs,"\n")
        # cat(loq_obs_master,"\n")
        # if(sum(loq_obs_master==2)==1) browser()
        # by rows!!
        lloq_mat = np.repeat(loq_full[loq_obs_master == 2], loq_obs.shape[0])
        lloq_mat = lloq_mat.reshape(loq_obs.shape[0], (lloq_mat.size)/loq_obs.shape[0])
        uloq_mat = np.repeat(uloq_full[loq_obs_master == 2], loq_obs.shape[0])
        uloq_mat = uloq_mat.reshape(loq_obs.shape[0], (uloq_mat.size)/loq_obs.shape[0])
        
        
        loq_comb_l = loq_obs*np.nan
        loq_comb_u = loq_obs*np.nan
        
        # BLOQ
        loq_comb_l[loq_obs == 1] = -np.Infinity
        loq_comb_u[loq_obs == 1] = lloq_mat[loq_obs == 1]
        
        # ULOQ
        loq_comb_l[loq_obs == 2] = uloq_mat[loq_obs == 2]
        loq_comb_u[loq_obs == 2] = np.Infinity
        
        # normal observations
        loq_comb_l[loq_obs == 0] = lloq_mat[loq_obs == 0]
        loq_comb_u[loq_obs == 0] = uloq_mat[loq_obs == 0]
        
        # compute all probabilities
        pred_pot_loq = pred[loq_obs_master == 2]
        cov_pot_loq = cov[loq_obs_master == 2,loq_obs_master == 2]
        p_loq_comb = np.repeat(0, loq_obs.shape[0])
        p_loq_comb_full = np.repeat(0, loq_obs.shape[0]) # for diagnostics
        for j in range(0,loq_obs.shape[0]):
            p_loq_comb_tmp = stats.mvn.mvnun(loq_comb_l[j,:], loq_comb_u[j,:], pred_pot_loq, cov_pot_loq)
            #p_bloq_comb_tmp = mnormt::sadmvn(bloq_comb_l[j,],bloq_comb_u[j,], pred, cov)
        
            # filter out low probability values
            p_loq_comb_full[j] = p_loq_comb_tmp # save initial probs for diagnostics
            if p_loq_comb_tmp < loq_prob_limit: 
                p_loq_comb_tmp = 0
            p_loq_comb[j] = p_loq_comb_tmp
        
        
        # sum of probabilities
        tot_p = sum(p_loq_comb_full)
        max_diff = PI_alpha/2*len(loq_obs_master == 2) # max p missed if all points are truncated with PI 
        if tot_p > 1.01 or tot_p < (1-max_diff):
            raise Exception("Sum of initial probabilities: %6.5g\n" + "Probabilities do not add up to one!" % tot_p)
        
        # rescale probabilities
        p_loq_comb = p_loq_comb/sum(p_loq_comb)
        
        if verbose is True:
            loq_obs_tmp = loq_obs_master 
            xt = None
            model = xt
            for j in range(0, loq_obs.shape[0]):
                loq_obs_tmp[loq_obs_master==2] = loq_obs[j,:] 
                df_p = pd.DataFrame({"model": [model_switch_i],"xt": [xt_i], "pred": [pred], "LOQ":[loq_obs_tmp]})
                df_p = df_p.sort_values(by=["model", "xt"])
                #print(df_p)
                print("Time: %1.f" + "\nLOQ: %1.f" + "\np_initial: %8.4g" + "\np_final: %8.4g" + "\n\n" % (df_p["xt"], df_p["LOQ"], p_loq_comb_full[j], p_loq_comb[j]))
            
            print("sum of initial probabilities: %6.5g\n" % tot_p)
            print("sum of final probabilities: %6.5g\n" % sum(p_loq_comb))
            print("\n")
        
        
        # compute FIM for each case and combine
        fim = zeros(fim_size)
        loq_obs_tmp = loq_obs_master 
        for j in range(0, loq_obs.shape[0]):
            #j=2
            loq_obs_tmp[loq_obs_master==2] = loq_obs[j,:] 
            if any(loq_obs_tmp==0) and p_loq_comb[j]!=0:
                fim_tmp = mf_all(model_switch_i[loq_obs_tmp==0,1], xt_i[loq_obs_tmp==0,1], x_i,a_i,bpop_val,d_full,sigma_full,docc_full,poped_db)["ret"]
                fim = fim + p_loq_comb[j]*fim_tmp


    else: # D2 method for BLOQ
        fim = zeros(fim_size)
        if any(loq_obs_master==0):
            fim = mf_all(model_switch_i[loq_obs_master==0,1], xt_i[loq_obs_master==0,1], x_i,a_i,bpop_val,d_full,sigma_full,docc_full,poped_db)["ret"]
                
    
    return {"fim": fim, "poped_db": poped_db}
