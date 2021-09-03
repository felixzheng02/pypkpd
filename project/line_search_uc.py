"""

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from project.util import trans
from project.feval import feval

def line_search_uc(x0, f0, g0, d, f_handle,f_options,exp_index,
                           ls_stepmax=1, #max step length for line search
                           ls_delta_alpha=1e-4, #convergence criterion on alpha
                           ls_fdecreas=1e-4, #sufficent decrease for line search
                           allow_negative=False
                           ):
    dim = x0.size
    if np.linalg.norm(d) > ls_stepmax:
        d = d*ls_stepmax/np.linalg.norm(d)
    
    #compute lambda_min
    lambda_min = ls_delta_alpha/max(abs(d)/max(np.array([[abs(x0)], [1,dim,1]]).reshape[1,2]))
    lambda_num = 1 #newton step
    slope = np.matmul(np.transpose(g0), d) #slope at beginning of line search
    
    #trans  = function(x) matrix(c(x[bpop_index],exp(x[d_index])),ncol=1,byrow=T)
    
    
    while True:
        # compute new x
        x1 = x0 + lambda_num*d
        
        # if x cannot be negative
        if allow_negative is False:
            for i in range(0, x1.shape[0]):
                for j in range(0, x1.shape[1]):
                    if x1[i,j] < 0: 
                        x1[i,j] = 0.00001
                    else:
                        x1[i,j] = x1[i,j]
                    j = j + 1
                i = i + 1
        
        #evaluate function
        f_options[[1]] = trans(x1)

        if any(np.isnan(f_options[[1]])): 
            browser()
        
        fval1 = feval(f_handle,f_options)
        fval1 = fval1["k"]
        if lambda_num < lambda_min:
            x_min = x0
            f_min = f0
            break
        elif fval1 < f0+ls_fdecreas*lambda_num*slope:
            x_min = x1
            f_min = fval1
            break
        else:
            if lambda_num ==1:
                lambda_tmp = -slope/(2*(fval1-f0-slope))
            else:        
                rhs1 = fval1-f0-lambda_num*slope
                rhs2 = fval2-f0-lambda_num2*slope
                a = (rhs1/lambda_num^2-rhs2/lambda_num2^2)/(lambda_num-lambda_num2)
                b = (-lambda_num2*rhs1/lambda_num^2+lambda_num*rhs2/lambda_num2^2)/(lambda_num-lambda_num2)
                if a == 0:
                    lambda_tmp=-slope/(2*b)
                else:
                    disc = b^2-3*a*slope
                    if disc < 0:
                        lambda_tmp = 0.5*lambda_num
                    elif b < 0:
                        lambda_tmp = (-b*np.sqrt(disc))/(3*a)
                    else:
                        lambda_tmp = -slope/(b+np.sqrt(disc))

                    lambda_tmp = min(np.array([0.5*lambda_num, lambda_tmp]).reshape[1,2])
                

        lambda_num2 = lambda_num
        lambda_num = max(lambda_tmp, 0.1*lambda_num)
        fval2 = fval1
    
    return {"x_min": x_min, "f_min": f_min}


