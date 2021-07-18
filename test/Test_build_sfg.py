



import numpy as np
from project.build_sfg import build_sfg
from project.models import ff_PK_1_comp_oral_md_CL



build_sfg(model=ff_PK_1_comp_oral_md_CL)
#> function (x, a, bpop, b, bocc) 
#> parameters = c(CL = bpop[1] * exp(b[1]), Favail = bpop[2], KA = bpop[3] * 
#>     exp(b[2]), V = bpop[4] * exp(b[3]), DOSE = a[1], TAU = a[2])
#> <environment: 0x7fe20d5952a8>
etas = np.array({"Favail": "exp", "KA": "exp", "V": "add", "CL": "exp"})
build_sfg(model=ff_PK_1_comp_oral_md_CL, etas=etas)
#> function (x, a, bpop, b, bocc) 
#> parameters = c(CL = bpop[1] * exp(b[1]), Favail = bpop[2] * exp(b[2]), 
#>     KA = bpop[3] * exp(b[3]), V = bpop[4] + b[4], DOSE = a[1], 
#>     TAU = a[2])
#> <environment: 0x7fe20d5952a8>