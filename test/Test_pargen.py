"""
############# START #################
## Create PopED database
## (warfarin example)
#####################################

## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## find the parameters that are needed to define from the structural model
##ff.PK.1.comp.oral.sd.CL
#> function (model_switch, xt, parameters, poped_db) 
#> {
#>     with(as.list(parameters), {
#>         y = xt
#>         y = (DOSE * Favail * KA/(V * (KA - CL/V))) * (exp(-CL/V * 
#>             xt) - exp(-KA * xt))
#>         return(list(y = y, poped_db = poped_db))
#>     })
#> }
#> <bytecode: 0x7fe20a979808>
#> <environment: namespace:PopED>
## -- parameter definition function 
## -- names match parameters in function ff
"""

import numpy as np
import scipy as sp
from project.ones import ones
from project.pargen import pargen
from project.models import feps_prop
from project.models import ff_PK_1_comp_oral_sd_CL
from project.create_poped_database import sfg
from project.create_poped_database import create_poped_database

## -- Define model, parameters, initial design
poped_db = create_poped_database(ff_fun=ff_PK_1_comp_oral_sd_CL,
                                  fg_fun=sfg,
                                  fError_fun=feps_prop,
                                  bpop=np.array({"CL": 0.15, "V": 8, "KA": 1.0, "Favail": 1}), 
                                  notfixed_bpop=np.array([1,1,1,0]),
                                  d=np.array({"CL": 0.07, "V": 0.02, "KA": 0.6}), 
                                  sigma={"prop": 0.01},
                                  groupsize=32,
                                  xt=np.array([0.5,1,2,6,24,36,72,120]),
                                  a={"DOSE": 70})

print(poped_db)
############# END ###################
## Create PopED database
## (warfarin example)
#####################################


# Adding 40% Uncertainty to fixed effects log-normal (not Favail)
bpop_vals = np.array({"CL": 0.15, "V": 8, "KA": 1.0, "Favail": 1})
bpop_vals_ed_ln = np.stack((ones(bpop_vals.size, 1)*4, # log-normal distribution
                      bpop_vals,
                      ones(bpop_vals.size, 1)*(bpop_vals*0.4)^2), axis=1) # 40% of bpop value
bpop_vals_ed_ln["Favail",:]  = np.array([0,1,0])

pars_ln = pargen(par=bpop_vals_ed_ln,
               user_dist_pointer=None,
               sample_size=1000,
               bLHS=1,
               sample_number=None,
               poped_db=poped_db)


# Adding 10% Uncertainty to fixed effects normal-distribution (not Favail)
bpop_vals_ed_n = np.stack((ones(bpop_vals.size, 1)*1, # log-normal distribution
                      bpop_vals,
                      ones(bpop_vals.size, 1)*(bpop_vals*0.1)^2), axis=1) # 10% of bpop value
bpop_vals_ed_n["Favail",:]  = np.array([0,1,0])

pars_n = pargen(par=bpop_vals_ed_n,
               user_dist_pointer=None,
               sample_size=1000,
               bLHS=1,
               sample_number=None,
               poped_db=poped_db)


# Adding 10% Uncertainty to fixed effects uniform-distribution (not Favail)
bpop_vals_ed_u = np.stack((ones(bpop_vals.size, 1)*2, # uniform distribution
                        bpop_vals,
                        ones(bpop_vals.size, 1)*(bpop_vals*0.1)), axis=1) # 10% of bpop value
bpop_vals_ed_u["Favail",:]  = np.array([0,1,0])

pars_u = pargen(par=bpop_vals_ed_u,
                 user_dist_pointer=None,
                 sample_size=1000,
                 bLHS=1,
                 sample_number=None,
                 poped_db=poped_db)


# Adding user defined distributions
bpop_vals_ed_ud = np.stack((ones(bpop_vals.size, 1)*3, # user dfined distribution
                         bpop_vals,
                         bpop_vals*0.1), axis=1) # 10% of bpop value
bpop_vals_ed_ud["Favail",:]  = np.array([0,1,0])

# A normal distribution
def my_dist(*args):
    par_vec = np.random.normal(bpop_vals_ed_ud[:,1], bpop_vals_ed_ud[:,2], np.array([1,1,1,1]))


pars_ud = pargen(par=bpop_vals_ed_ud,
                  user_dist_pointer=my_dist,
                  sample_size=1000,
                  bLHS=1,
                  sample_number=None,
                  poped_db=poped_db)


