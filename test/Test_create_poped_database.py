## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

import numpy as np
#from numpy.core.records import array
#from project.models import feps_prop
# from project.evaluate_design import evaluate_design

from project.models import feps_prop
from project.models import ff_PK_1_comp_oral_sd_CL
from project.create_poped_database import sfg
from project.create_poped_database import create_poped_database


## -- Define initial design  and design space
poped_db = create_poped_database(
                                 ff_file=ff_PK_1_comp_oral_sd_CL,
                                 fg_file=sfg,
                                 fError_file=feps_prop,
                                 bpop=np.array([0.15, 8, 1.0, 1]), 
                                 notfixed_bpop=np.array([1,1,1,0]),
                                 d=np.array([0.07, 0.02, 0.6]), 
                                 sigma=0.01,
                                 groupsize=32,
                                 xt=np.array([0.5,1,2,6,24,36,72,120]),
                                 minxt=0,
                                 maxxt=120,
                                 a=70)
print(poped_db)
#> Warning: cannot open file 'sfg': No such file or directory#> Error in file(filename, "r", encoding = encoding): cannot open the connection

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



## evaluate initial design
# evaluate_design(poped_db)
#> Error in calc_ofv_and_fim(poped_db, ...): object 'poped_db' not found