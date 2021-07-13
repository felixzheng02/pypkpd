"""
## Get size of the fisher information matrix


## Author: Caiya Zhang, Yuchen Zheng
"""




def get_fim_size(poped_db):
    #Returns the size of FIM, i$e. col or row size
    numnotfixed_bpop = sum(poped_db["parameters"]["notfixed_bpop"])
    numnotfixed_d    = sum(poped_db["parameters"]["notfixed_d"])
    numnotfixed_covd = sum(poped_db["parameters"]["notfixed_covd"])
    numnotfixed_docc  = sum(poped_db["parameters"]["notfixed_docc"])
    numnotfixed_covdocc  = sum(poped_db["parameters"]["notfixed_covdocc"])
    numnotfixed_sigma  = sum(poped_db["parameters"]["notfixed_sigma"])
    numnotfixed_covsigma  = sum(poped_db["parameters"]["notfixed_covsigma"])
    
    n = numnotfixed_bpop + numnotfixed_d + numnotfixed_covd + numnotfixed_docc + numnotfixed_covdocc + numnotfixed_sigma + numnotfixed_covsigma
    
    return n

