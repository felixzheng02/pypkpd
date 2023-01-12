import path
from project.model_prediction import model_prediction
from project.sfg import Sfg
from project.models import Models
from project.create_poped_database import create_poped_database


ff_PK_1_comp_oral_sd_CL = Models.ff_PK_1_comp_oral_sd_CL
feps_prop = Models.feps_prop
sfg = Sfg.sfg_2

pypkpd_db = create_poped_database(
                                    ff_fun = ff_PK_1_comp_oral_sd_CL,
                                    fg_fun = sfg,
                                    fError_fun = feps_prop
                                    )