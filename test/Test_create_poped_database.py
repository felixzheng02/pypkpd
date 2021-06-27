## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

import numpy as np
from numpy.core.records import array
from project.evaluate_design import evaluate_design
from project.create_poped_database import create_poped_database

## find the parameters that are needed to define from the structural model
## ff.PK.1.comp.oral.md.CL

#> function (model_switch, xt, parameters, poped_db) 
#> {
#>     with(as.list(parameters), {
#>         y = xt
#>         N = floor(xt/TAU) + 1
#>         y = (DOSE * Favail/V) * (KA/(KA - CL/V)) * (exp(-CL/V * 
#>             (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - 
#>             exp(-CL/V * TAU)) - exp(-KA * (xt - (N - 1) * TAU)) * 
#>             (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
#>         return(list(y = y, poped_db = poped_db))
#>     })
#> }
#> <bytecode: 0x7fe25cadfcb0>
#> <environment: namespace:PopED>
## -- parameter definition function 
## -- names match parameters in function ff
def sfg (x,a,bpop,b,bocc):
    paras_dict = {"CL":bpop[1]*np.exp(b[1]),
                "V":bpop[2]*np.exp(b[2]),
                "KA":bpop[3]*np.exp(b[3]),
                "Favail":bpop[4],
                "DOSE":a[1]}
    parameters = np.array([paras_dict])
    return parameters


## -- Define initial design  and design space
poped_db = create_poped_database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                  fg_file="sfg",
                                  fError_file="feps.prop",
                                  bpop=np.array([0.15, 8, 1.0, 1]), 
                                  notfixed_bpop=np.array([1,1,1,0]),
                                  d=np.array([0.07, 0.02, 0.6]), 
                                  sigma=0.01,
                                  groupsize=32,
                                  xt=np.array([0.5,1,2,6,24,36,72,120]),
                                  minxt=0,
                                  maxxt=120,
                                  a=70)
#> Warning: cannot open file 'sfg': No such file or directory#> Error in file(filename, "r", encoding = encoding): cannot open the connection

## evaluate initial design
evaluate_design(poped_db)
#> Error in calc_ofv_and_fim(poped_db, ...): object 'poped_db' not found