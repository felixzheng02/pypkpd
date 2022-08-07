"""
## Evaluate a design
## 
## This function evaluates the design defined in a poped database.
## 
## @param pypkpd_db A poped database
## @param ... Extra parameters passed to \code{\link{calc_ofv_and_fim}} and \code{\link{get_rse}}
## @return A list of elements evaluating the current design.
## @export
## 
## @family evaluate_design

## Author: Caiya Zhang, Yuchen Zheng
"""



from project.get_cv import get_rse
from project.calc_ofv_and_fim import calc_ofv_and_fim
from project.util import get_dict_value

def evaluate_design(pypkpd_db, **kwargs):
    out = calc_ofv_and_fim(pypkpd_db, **kwargs)
    if get_dict_value(out, "fim") is None:
        out["rse"] = None
    else:
        out["rse"] = get_rse(get_dict_value(out, "fim"), pypkpd_db, **kwargs)
    
    # out["fim"] datanam should be colnam and rownam
    
    return out