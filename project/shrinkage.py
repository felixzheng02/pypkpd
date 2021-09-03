"""
#' Predict shrinkage of empirical Bayes estimates (EBEs) in a population model
#'
#' @param poped_db database
#' @param num_sim_ids If \code{use_mc=TRUE}, how many individuals should be
#'   simulated to make the computations.
#' @param use_mc Should the calculation be based on monte-carlo simulations. If
#'   not then then a first order approximation is used
#' @param use_purrr If \code{use_mc=TRUE} then should the method use the package
#'   purrr in calculations?  This may speed up computations (potentially).
#'
#' @return The shrinkage computed in variance units, standard deviation units
#'   and the relative standard errors of the EBEs.
#' @export
#'
#' @references \enumerate{ 
#'   \item Combes, F. P., Retout, S.,
#'   Frey, N., & Mentre, F. (2013). Prediction of shrinkage of individual
#'   parameters using the Bayesian information matrix in non-linear mixed effect
#'   models with evaluation in pharmacokinetics. Pharmaceutical Research, 30(9),
#'   2355-67. \doi{10.1007/s11095-013-1079-3}. 
#'   \item Hennig, S., Nyberg, J., Fanta, S., Backman, J.
#'   T., Hoppu, K., Hooker, A. C., & Karlsson, M. O. (2012). Application of the
#'   optimal design approach to improve a pretransplant drug dose finding design
#'   for ciclosporin. Journal of Clinical Pharmacology, 52(3), 347-360.
#'   \doi{10.1177/0091270010397731}. 
#'   }

## Author: Caiya Zhang, Yuchen Zheng
"""

import re
import inspect
import warnings
import numpy as np
from project.size import size
from project.zeros import zeros
from matpy.matrix import matrix
from project.get_cv import get_rse
from project.getfulld import getfulld
from project.get_parnam import get_parnam
from project.evaluate_fim import evaluate_fim
from project.find_largest_index import find_largest_index
from project.create_poped_database import create_poped_database


def shrinkage(poped_db,
                use_mc = False,
                num_sim_ids = 1000,
                use_purrr = False):
  
    # if (poped_db["design"]m > 1) {
    #   warning("Shrinkage should only be computed for a single arm, please adjust your script accordingly.")
    # }
    
    group = typ = None
    
    # tranform random effects to fixed effects 
    tmp_fg = poped_db["model"]["fg_pointer"]
    if type(tmp_fg) is str:
        tmp_fg = eval(tmp_fg)
    largest_bpop = find_largest_index(tmp_fg,lab = "bpop")
    largest_b = find_largest_index(tmp_fg,lab = "b")
    largest_eps = find_largest_index(poped_db["model"]["ferror_pointer"], lab="epsi", mat=True, mat_row=False)
    def replace_fun(txt):
        #extract number
        num = float(re.findall(txt, "\d+")) %>% "+"(.,largest_bpop)
        txt = txt.replace("\\d+", num)
    
    txt = inspect.getsource(body(tmp_fg))
    txt = txt.replace("b\\[(\d+)", replace_fun)
    txt = txt.replace("b\\[", "bpop\\[")
    
    env = environment()
    body(tmp_fg, envir=env) = parse(text=txt)
    #environment(sfg_tmp) = env
    
    # fix population parameters
    # add IIV and IOV as fixed effects
    # change to one individual
    # TODO: need to add iov
    extra_bpop = matrix(0, (largest_b, 3))
    bpop_tmp = matrix([poped_db["parameters"]["bpop"], extra_bpop])
    notfixed_d = poped_db["parameters"]["notfixed_d"]
    non_zero_d = poped_db["parameters"]["d"][:,1] != 0
    
    notfixed_d_new = notfixed_d or non_zero_d
    
    notfixed_bpop_tmp = matrix([np.repeat(0, largest_bpop), notfixed_d_new])
    poped_db_sh = create_poped_database(poped_db,
                                        fg_fun = tmp_fg,
                                        bpop = bpop_tmp, 
                                        nbpop = size(bpop_tmp)[0],
                                        notfixed_bpop = notfixed_bpop_tmp,
                                        d = matrix(nrow=0, ncol=3),
                                        notfixed_d = matrix(nrow = 0, ncol = 3),
                                        NumRanEff = 0,
                                        notfixed_sigma = matrix(np.repeat(0, largest_eps)),
                                        groupsize = 1,
                                        mingroupsize = 1,
                                        mintotgroupsize = 1)
    
    # Compute FIM_map 
    fulld = getfulld(poped_db["parameters"]["d"][notfixed_d_new,1], poped_db["parameters"]["covd"])
    
    # get just groupwise values as well
    db_list = list()
    #db_list[["all"]] = poped_db_sh
    num_groups = poped_db_sh["design"]["m"]
    for i in range(0, num_groups):
        tmp_design = poped_db_sh["design"]
        tmp_design["m"] = 1
        for nam in tmp_design.get_datanam()[tmp_design.get_datanam() != "m"]:
            tmp_design[nam] = tmp_design[nam][i,:]
        

        tmp_db = poped_db_sh
        tmp_db["design"] = tmp_design
        
        db_list["grp_"+i] = tmp_db
    
    
    #out_list = list()
    out_df = matrix()
    #for(i in 1:1){
    for i in range(0, len(db_list)):
        
        poped_db_sh = db_list[i]
        
        if use_mc is True:
            # create list of eta values from pop model 
            if any(size(fulld) ==0):
                b_sim_matrix = zeros(num_sim_ids,0)
            else:
                b_sim_matrix = np.random.multivariate_normal(num_sim_ids, sigma=fulld)            
            
            
            # create list of occasion eta values from pop model 
            NumOcc=poped_db["parameters"]["NumOcc"]
            if NumOcc != 0: 
                warnings.warn("Shrinkage not computed for occasions\n")
            fulldocc = getfulld(poped_db["parameters"]["docc"][:,1], poped_db["parameters"]["covdocc"])
            bocc_sim_matrix = zeros(num_sim_ids*NumOcc, len(poped_db["parameters"]["docc"][:,1]))
            if size(fulldocc)[0] != 0: 
                bocc_sim_matrix = np.random.multivariate_normal(num_sim_ids*NumOcc, sigma=fulldocc)
            
            #now use these values to compute FIMmap
            if use_purrr is True:
                fim_mean = matrix(np.transpose(b_sim_matrix)) %>% as.list() %>% 
                purrr::map(function(x){
                    x_tmp = matrix(0,nrow = largest_b,ncol = 3)
                    x_tmp[,2] = t(x)
                    bpop_tmp =rbind(poped_db["parameters"]bpop,x_tmp)
                    poped_db_sh_tmp = create.poped.database(poped_db_sh, bpop=bpop_tmp)
                    evaluate.fim(poped_db_sh_tmp)
                }) %>% simplify2array() %>% apply(1:2, mean)
            else:
                fim_tmp = 0
                for i in range(0, num_sim_ids):
                    x_tmp = matrix(0, (largest_b, 3))
                    x_tmp[:,1] = np.transpose(b_sim_matrix[i,:])
                    bpop_tmp = matrix([poped_db["parameters"]["bpop"], x_tmp])
                    poped_db_sh_tmp = create_poped_database(poped_db_sh, bpop=bpop_tmp)
                    fim_tmp = fim_tmp + evaluate_fim(poped_db_sh_tmp)
                
                fim_mean = fim_tmp/num_sim_ids
             
            fim_map = fim_mean + np.linalg.inv(fulld)
        else:
            fim_map = evaluate_fim(poped_db_sh) + np.linalg.inv(fulld)
        
        
        rse = get_rse(fim_map,poped_db = poped_db_sh)
        poped_db_tmp = poped_db
        poped_db_tmp["parameters"]["notfixed_d"] = notfixed_d_new
        parnam = get_parnam(poped_db_tmp)
        names(rse) = parnam[re.findall("^(D\\[|d\\_)", parnam)]
        
        #names(rse) = paste0("d[",1:length(rse),"]")
        
        # shrinkage on variance scale
        shrink = np.diag(np.matmul(np.linallg.inv(fim_map), np.linag.inv(fulld)))
        names(shrink) =  names(rse)
        # shrinkage on SD scale
        var_n = (1-shrink)*np.diag(fulld)
        shrink_sd = 1-np.sqrt(var_n)/np.sqrt(np.diag(fulld))
        names(shrink_sd) =  names(rse)
        

        out_df_tmp = matrix([shrink, shrink_sd, rse])
        out_df_tmp["type"] = matrix(["shrink_var","shrink_sd","se"])
        out_df_tmp["group"] = matrix(db_list[i].get_datanam())
        out_df = matrix([out_df, out_df_tmp])
        #out_list[[names(db_list[i])]] = list(shrinkage_var=shrink,shrinkage_sd=shrink_sd, rse=rse)
    
    
    #return(list(shrinkage_var=shrink,shrinkage_sd=shrink_sd, rse=rse))
    #return(out_list)
    out_df = np.arrange(out_df, dplyr::desc(type), group)
    
    if poped_db["design"]["m"] > 1:
        weights = poped_db["design"]["groupsize"]/sum(poped_db["design"]["groupsize"])
        data_tmp = out_df 
        data_tmp[data_tmp == 1] = np.nan 
        
        data_tmp = data_tmp %>% dplyr::group_by(type) %>% 
            dplyr::summarise_at(dplyr::vars(dplyr::starts_with(c('d[','d_'))),
                                list(~ stats::weighted.mean(., weights,na.rm = T)))
            
        data_tmp["group"] = "all_groups"
        out_df = matrix([out_df, data_tmp])
        out_df = dplyr::arrange(out_df,dplyr::desc(type),group)
    
    return out_df
  

