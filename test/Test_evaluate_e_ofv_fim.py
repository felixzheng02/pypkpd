"""
## library(PopED)

############# START #################
## Create PopED database
## (warfarin model for optimization
##  with parameter uncertainty)
#####################################

## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error
## to avoid sample times at very low concentrations (time 0 or very late samoples).

## find the parameters that are needed to define from the structural model


## Author: Caiya Zhang, Yuchen Zheng
"""

import numpy as np
from numpy.core.records import array
from project.tic import tic
from project.toc import toc
from project.ones import ones
from project.evaluate_e_ofv_fim import evaluate_e_ofv_fim
from project.create_poped_database import create_poped_database



def sfg(x,a,bpop,b,bocc):
    CL = bpop[1]*np.exp(b[1]),
    V = bpop[2]*np.exp(b[2]),
    KA = bpop[3]*np.exp(b[3]),
    Favail = bpop[4],
    DOSE = a[1]
    parameters = np.array([CL,V,KA,Favail,DOSE])
    return parameters


# Adding 10% log-normal Uncertainty to fixed effects (not Favail)
CL = 0.15
V = 8
KA = 1.0
Favail = 1
bpop_vals = np.array([CL,V,KA,Favail])
bpop_vals_ed_ln = np.concatenate(ones(bpop_vals.size,1)*4, bpop_vals,ones(bpop_vals.size,1)*(bpop_vals*0.1)^2) ## log-normal distribution, 10% of bpop value
bpop_vals_ed_ln["Favail",:]  = np.array(0,1,0)
bpop_vals_ed_ln
#>          bpop_vals         
#> CL     4      0.15 0.000225
#> V      4      8.00 0.640000
#> KA     4      1.00 0.010000
#> Favail 0      1.00 0.000000
## -- Define initial design  and design space
poped_db = create_poped_database(ff_fun="ff_PK_1_comp_oral_sd_CL",
                                fg_fun=sfg(),
                                fError_fun="feps_add_prop",
                                bpop=bpop_vals_ed_ln, 
                                notfixed_bpop=np.array([1,1,1,0]),
                                #CL=0.07, V=0.02, KA=0.6
                                d=np.array([0.07, 0.02, 0.6]), 
                                sigma=np.array([0.01,0.25]),
                                groupsize=32,
                                xt=np.array([0.5,1,2,6,24,36,72,120]),
                                minxt=0,
                                maxxt=120,
                                a=70,
                                mina=0,
                                maxa=100)

############# END ###################
## Create PopED database
## (warfarin model for optimization
##  with parameter uncertainty)
#####################################


## ED evaluate (with very few samples)
output = evaluate_e_ofv_fim(poped_db,ED_samp_size=10)
output["E_ofv"]
#> [1] 55.45214
## API evaluate (with very few samples)
output = evaluate_e_ofv_fim(poped_db,ED_samp_size=10,ofv_calc_type=4)
output["E_ofv"]
#> [1] 55.46088
## ED evaluate using Laplace approximation 
tic()
output = evaluate_e_ofv_fim(poped_db,use_laplace=True)
toc()
#> Elapsed time: 1.3 seconds.output$E_ofv
#> [1] 1.302806e+24
if False:

    ## ED expected value with more precision. 
    ## Compare time and value to Laplace approximation.
    ## Run a couple of times to see stochasticity of calculation.
    tic()
    e_ofv_mc = evaluate_e_ofv_fim(poped_db,ED_samp_size=500)
    toc()
    e_ofv_mc["E_ofv"]
    
    # If you want to get an E(FIM) from the laplace approximation you have to ask for it
    # and it will take more time.
    output = evaluate_e_ofv_fim(poped_db,use_laplace=True,laplace_fim=True)
    output["E_fim"]
  
 


