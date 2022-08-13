"""
## Parameter simulation
## 
## Function generates random samples for a list of parameters
## Function similar to that in matlab
## 
## @param par A Matrix describing the parameters. Each row is a parameter and 
##   the Matrix has three columns: 
##   \enumerate{ 
##   \item First column - Type of
##   distribution (0-fixed, 1-normal, 2-uniform, 3-user specified, 4-lognormal,
##   5-Truncated normal). 
##   \item Second column - Mean of distribution. 
##   \item Third
##   column - Variance or range of distribution. 
##   }
## @param user_dist_pointer A text string of the name of a function that
##   generates random samples from a user defined distribution.
## @param sample_size The number of random samples per parameter to generate
## @param bLHS Logical, indicating if Latin Hypercube Sampling should be used.
## @param sample_number The sample number to extract from a user distribution.
## @param pypkpd_db A PopED database.
##   
## @return A Matrix of random samples of size (sample_size x number_of_parameters)
## @example test/Test_pargen.py
## @export


## Author: Caiya Zhang, Yuchen Zheng
"""


from scipy import stats
from scipy.stats import norm
import numpy as np
# from scipy.linalg import norm
from project.size import size
from matpy.matrix import Matrix
from project.feval import feval
from project.getTruncatedNormal import getTruncatedNormal


def pargen (par: Matrix,
            user_dist_pointer,
            sample_size,
            bLHS,
            sample_number,
            pypkpd_db):
    
    nvar = par.get_shape(1)
    ret = Matrix(np.zeros(sample_size, nvar))
    
    # for log-normal distributions
    # mu=log(par[,2]^2/sqrt(par[,3]+par[,2]^2))
    # sd=sqrt(log(par[,3]/par[,2]^2+1))
    # exp(rnorm(100000,mu,si))
    # mean(exp(rnorm(100000,mu,si)))
    # exp(mu+qnorm(P)*sd)) 
    # exp((log(par[,2]^2/sqrt(par[,3]+par[,2]^2)))+qnorm(P)*(sqrt(log(par[,3]/par[,2]^2+1)))) 
    
    ## using rlnorm (may be faster)
    # Adding 40% Uncertainty to fixed effects log-normal (not Favail)
    # bpop_vals = c(CL=0.15, V=8, KA=1.0, Favail=1)
    # bpop_vals_ed_ln = cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
    #                          bpop_vals,
    #                          ones(length(bpop_vals),1)*(bpop_vals*0.4)^2) # 40% of bpop value
    # bpop_vals_ed_ln["Favail",]  = c(0,1,0)
    # 
    # # with log-normal distributions
    # pars.ln = pargen(par=bpop_vals_ed_ln,
    #                   user_dist_pointer=None,
    #                   sample_size=1000,
    #                   bLHS=1,
    #                   sample_number=None,
    #                   pypkpd_db)
    # sample_size=1000
    # data = apply(bpop_vals_ed_ln,1,
    #               function(par) rlnorm(sample_size,log(par[2]^2/sqrt(par[3]+par[2]^2)),
    #                                    sqrt(log(par[3]/par[2]^2+1))))
    # colMeans(data)
    # var(data)
    
    
    if bLHS == 0: #Random Sampling
        for k in range(0, sample_size):
            np = par.get_shape(0)
            if np != 0:
                n = np.random.randn(np, 1) # normal mean=0 and sd=1
                u = np.random.rand(np, 1) * 2 - 1 # uniform from -1 to 1
                t = par.get_partial_matrix([[None, None], [0, 1]]) # type of distribution
                c2 = par.get_partial_matrix([[None, None], [2, 3]]) # variance or range of distribution
                
                bUserSpecifiedDistribution: bool = sum(t.get_data() == 3) >= 1 #If at least one user specified distribution, sum(t == 3) counts number of 3 in t
                ret.data[k, :] = (t.get_data() == 0) * par.get_data()[:, 1] + (t.get_data() == 2) * (par.get_data()[:, 1] 
                + u * c2.get_data()/2) + (t.get_data() == 1) * (par.get_data()[:, 1] + n * c2.get_data()^(1/2))
                + (t.get_data() == 4) * np.exp((np.log(par.get_data()[:, 1]^2 / np.sqrt(par.get_data()[:, 2] 
                + par.get_data()[:, 1]^2))) + n * (np.sqrt(np.log(par.get_data()[:, 2]/par.get_data()[:, 1]^2+1))))
                #(t==4)*par[,2,drop=F]*exp(n*c2^(1/2))
                if sum(t == 5) > 0: #Truncated normal
                    for i in range(0, par.get_shape(0)):
                        if t.get_one_data([i, 0]) == 5:
                            ret.data[k, i] = getTruncatedNormal(par.get_data()[i, 2], c2.get_data()[i, 0])
                
                if bUserSpecifiedDistribution:
                    if len(sample_number) == 0:
                        ret.data[k, :] = eval(str(user_dist_pointer) + "(" + str(ret[k, :],t,k,pypkpd_db) + ")")
                    else:
                        ret.data[k, :] = eval(str(user_dist_pointer) + "(" + str(ret[k, :],t,sample_number,pypkpd_db) + ")")

    elif nvar != 0: #LHS
        ran = np.random.rand(sample_size, nvar)
        # method of Stein
        for j in range(0, nvar):
            idx = np.random.random_sample(sample_size)
            P = (idx - ran[:, j])/sample_size 
                  # probability of the cdf
            returnArgs_list = [par.get_one_data([j, 1]),  #point
                                par.get_one_data([j, 1]) + norm.ppf(P)*np.sqrt(par.get_one_data([j, 3])), # normal
                                par.get_one_data([j, 1]) - par.get_one_data([j, 3])/2 + P*par.get_one_data([j, 3]), #uniform
                                ret.get_data()[:, j], #Do nothing
                                np.exp((np.log(par.get_one_data([j, 2])^2/np.sqrt(par.get_one_data([j, 3])+par.get_one_data([j, 2])^2)))+norm.ppf(P)*(np.sqrt(np.log(par.get_one_data([j, 3])/par.get_one_data([j, 2])^2+1))))] #log-normal 
                                #par[j,2]*exp(qnorm(P)*sqrt(par[j,3])) #log-normal
            returnArgs = returnArgs_list[par.get_one_data([j, 0]) + 1]
            
        
        if returnArgs is None:
            raise Exception("Unknown distribution for the inverse probability function used in Latin Hypercube Sampling")
        
        ret.data[:, j] = returnArgs
        
        bUserSpecifiedDistribution = sum(par.get_data()[:, 0] == 3) >= 1 #If at least one user specified distribution
        
        if bUserSpecifiedDistribution is True:
            for k in range(0, sample_size):
                if len(sample_number) == 0:
                    ret.data[k, :] = eval(str(user_dist_pointer) + "(" + str(ret.get_data()[k, :],par.get_data()[:, 0],k,pypkpd_db) + ")")
                else:
                    ret.data[k, :] = eval(str(user_dist_pointer) + "(" + str(ret.get_data()[k, :],par.get_data()[:, 0],sample_number,pypkpd_db) + ")")
    
    #return type: Matrix
    return ret