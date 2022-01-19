"""
## Evaluate a design
## 
## This function evaluates the design defined in a poped database.
## 
## @param poped_db A poped database
## @param ... Extra parameters passed to \code{\link{calc_ofv_and_fim}} and \code{\link{get_rse}}
## @return A list of elements evaluating the current design.
## @export
## 
## @family evaluate_design

## Author: Caiya Zhang, Yuchen Zheng
"""



from project.get_cv import get_rse
from project.calc_ofv_and_fim import calc_ofv_and_fim

def evaluate_design(poped_db, *argv):
    out = calc_ofv_and_fim(poped_db, *argv)
    if out["fim"] is None:
        out["rse"] = None
    else:
        out["rse"] = get_rse(out["fim"], poped_db, *argv)
    
    if out["rse"].keys() is not None:
        rownames(out["fim"]) = out["rse"].keys()
        colnames(out["fim"]) = out["rse"].keys()
    
    return out