"""
## library(PopED)

############# START #################
## Create PopED database
## (warfarin model for optimization)
#####################################

## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error  
## to avoid sample times at very low concentrations (time 0 or very late samples).

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.sd.CL
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

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.evaluate_fim import evaluate_fim
from project.create_poped_database import create_poped_database
from project.models import feps_add_prop
from project.models import ff_PK_1_comp_oral_sd_CL

def sfg(x,a,bpop,b,bocc):
    parameters=np.arrayCL=bpop[1]*np.exp(b[1]),
    V=bpop[2]*np.exp(b[2]),
    KA=bpop[3]*np.exp(b[3]),
    Favail=bpop[4],
    DOSE=a[1]
    return parameters


## -- Define initial design  and design space
poped_db = create_poped_database(ff_fun=ff_PK_1_comp_oral_sd_CL,
                                  fg_fun=sfg(),
                                  fError_fun=feps_add_prop,
                                  bpop=np.array([0.15,8,1.0,1]), 
                                  notfixed_bpop=np.array([1,1,1,0]),
                                  d=np.array([0.07,0.02,0.6]), 
                                  sigma=np.array([0.01,0.25]),
                                  groupsize=32,
                                  xt=np.array([0.5,1,2,6,24,36,72,120]),
                                  minxt=0.01,
                                  maxxt=120,
                                  a=np.array([70]),
                                  mina=np.array([0.01]),
                                  maxa=np.array([100]))

############# END ###################
## Create PopED database
## (warfarin model for optimization)
#####################################
print(poped_db)

## evaluate initial design 
FIM = evaluate_fim(poped_db) 
print(FIM)
#>             [,1]      [,2]      [,3]         [,4]         [,5]        [,6]
#> [1,] 17141.83891 20.838375 10.011000 0.000000e+00     0.000000  0.00000000
#> [2,]    20.83837 17.268051 -3.423641 0.000000e+00     0.000000  0.00000000
#> [3,]    10.01100 -3.423641 49.864697 0.000000e+00     0.000000  0.00000000
#> [4,]     0.00000  0.000000  0.000000 2.324341e+03     9.770352  0.03523364
#> [5,]     0.00000  0.000000  0.000000 9.770352e+00 19083.877564 11.72131703
#> [6,]     0.00000  0.000000  0.000000 3.523364e-02    11.721317 38.85137516
#> [7,]     0.00000  0.000000  0.000000 7.268410e+02  9656.158553 64.78095548
#> [8,]     0.00000  0.000000  0.000000 9.062739e+01   266.487127  2.94728469
#>              [,7]        [,8]
#> [1,]      0.00000    0.000000
#> [2,]      0.00000    0.000000
#> [3,]      0.00000    0.000000
#> [4,]    726.84097   90.627386
#> [5,]   9656.15855  266.487127
#> [6,]     64.78096    2.947285
#> [7,] 192840.20092 6659.569867
#> [8,]   6659.56987  475.500111get_rse(FIM,poped_db)
#>        CL         V        KA      d_CL       d_V      d_KA  sig_prop   sig_add 
#>  5.096246  3.031164 14.260384 29.761226 36.681388 26.748640 32.011719 25.637971 
print(np.linalg.det(FIM))
#> [1] 1.143859e+24ofv_fim(FIM,poped_db,ofv_calc_type=1) # det(FIM)
#> [1] 1.143859e+24ofv_fim(FIM,poped_db,ofv_calc_type=2) # 1/trace_matrix(inv(FIM))
#> [1] 9.127328ofv_fim(FIM,poped_db,ofv_calc_type=4) # log(det(FIM)) 
#> [1] 55.39645ofv_fim(FIM,poped_db,ofv_calc_type=6) # Ds with fixed effects as "important"
#> [1] 16.49204ofv_fim(FIM,poped_db,ofv_calc_type=6,
ds_index = np.array([1,1,1,0,0,0,1,1]) # Ds with random effects as "important"
#> [1] 21.23143ofv_fim(FIM,poped_db,ofv_calc_type=7) # 1/sum(get_rse(FIM,poped_db,use_percent=FALSE))
#> [1] 0.5772714
print(ds_index)