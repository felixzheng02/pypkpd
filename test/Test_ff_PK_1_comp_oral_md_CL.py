
import path
import numpy as np
from matpy.matrix import Matrix
from project.get_cv import get_rse
from project.evaluate_fim import evaluate_fim
from project.create_poped_database import sfg
from project.create_poped_database import create_poped_database
from project.plot_model_prediction import plot_model_prediction
from project.models import feps_add_prop
from project.models import ff_PK_1_comp_oral_md_CL


## -- parameter definition function 
## -- names match parameters in function ff
def sfg(x,a,bpop,b,bocc):
    parameters = Matrix(np.array([bpop.get_data()[0]*np.exp(b.get_data()[0]), 
                                bpop.get_data()[1]*np.exp(b.get_data()[1]),
                                bpop.get_data()[2]*np.exp(b.get_data()[2]),
                                bpop.get_data()[3],
                                a.get_data()[0],
                                a.get_data()[1]]),
                        (1, 6),
                        ["V", "KA", "CL", "Favail", "DOSE", "TAU"],
                        None, None)
    return parameters

## -- Define design and design space
poped_db = create_poped_database(ff_fun=ff_PK_1_comp_oral_md_CL,
                                  fg_fun=sfg,
                                  fError_fun=feps_add_prop,
                                  groupsize=20,
                                  m=2,
                                  sigma=Matrix(np.array([0.04,5e-6])),
                                  bpop=Matrix(np.array([72.8,0.25,3.75,0.9])), 
                                  d=Matrix(np.array([0.09,0.09,0.25^2])), 
                                  notfixed_bpop=Matrix(np.array([1,1,1,0])),
                                  notfixed_sigma=Matrix(np.array([0,0])),
                                  xt=Matrix(np.array([1,2,8,240,245])),
                                  minxt=Matrix(np.array([0,0,0,240,240])),
                                  maxxt=Matrix(np.array([10,10,10,248,248])),
                                  a=Matrix(np.array([20,40],[24,24])),
                                  bUseGrouped_xt=1,
                                  maxa=Matrix(np.array([200,24])),
                                  mina=Matrix(np.array([0,24])))

##  create plot of model without variability 
plot_model_prediction(poped_db)

## evaluate initial design
FIM = evaluate_fim(poped_db) 
print(FIM)

det_FIM = np.linalg.det(FIM)
print(det_FIM)

rse_FIM = get_rse(FIM, poped_db)
print(rse_FIM)
