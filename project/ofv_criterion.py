"""
## Normalize an objective function by the size of the FIM matrix
## 
## Compute a normalized OFV based on the size of the FIM matrix.  This value can then be used in 
## efficiency calculations. This is NOT the OFV used in optimization, see \code{\link{ofv_fim}}. 
## 
## @param ofv_f An objective function
## @param num_parameters The number of parameters to use for normalization
## @param poped_db a poped database
## @inheritParams ofv_fim
## 
## @return The specified criterion value.
## 
## @family FIM
##
## @example test/text_ofv_criterion.py
## 
## @export
# 
## Author: Caiya Zhang, Yuchen Zheng
"""

import sys
import numpy as np


def ofv_criterion(ofv_f, num_parameters, poped_db):
  
    #Input: the ofv
    #Return the single value that should be maximized

    ofv_calc_type=poped_db["settings"]["ofv_calc_type"]
    criterion_value = 0
    
    if ofv_calc_type == 0:
        criterion_value = ofv_f
    
    if ofv_calc_type == 1:  #D-Optimal Design
        criterion_value = ofv_f^(1/num_parameters)

    if ofv_calc_type == 4:  #D-Optimal Design
        criterion_value = np.exp(ofv_f)^(1/num_parameters)
    
    if ofv_calc_type == 2:  #A-Optimal Design
        criterion_value = ofv_f/num_parameters
        
    if ofv_calc_type == 3:  #S-Optimal Design
        try:
            print('Criterion for S-optimal design not implemented yet')
        except:
            sys.exit(1)
    
    if ofv_calc_type == 6:  #Ds-Optimal design
        criterion_value = ofv_f^(1/sum(poped_db["parameters"]["ds_index"]))
    
    
    return criterion_value
