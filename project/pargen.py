"""
#' Parameter simulation
#' 
#' Function generates random samples for a list of parameters
#' 
#' @param par A matrix describing the parameters. Each row is a parameter and 
#'   the matrix has three columns: 
#'   \enumerate{ 
#'   \item First column - Type of
#'   distribution (0-fixed, 1-normal, 2-uniform, 3-user specified, 4-lognormal,
#'   5-Truncated normal). 
#'   \item Second column - Mean of distribution. 
#'   \item Third
#'   column - Variance or range of distribution. 
#'   }
#' @param user_dist_pointer A text string of the name of a function that
#'   generates random samples from a user defined distribution.
#' @param sample_size The number of random samples per parameter to generate
#' @param bLHS Logical, indicating if Latin Hypercube Sampling should be used.
#' @param sample_number The sample number to extract from a user distribution.
#' @param poped_db A PopED database.
#'   
#' @return A matrix of random samples of size (sample_size x
#'   number_of_parameters)
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_pargen.R
#' @export
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work


## Author: Caiya Zhang, Yuchen Zheng
"""


from scipy import stats
from scipy.stats import norm
import numpy
# from scipy.linalg import norm
from project.size import size
from project.zeros import zeros
from project.feval import feval
from project.getTruncatedNormal import getTruncatedNormal


def pargen (par,user_dist_pointer,sample_size,bLHS,sample_number,poped_db):
  
    nvar = size(par)[0]
    ret = zeros(sample_size,nvar)
    
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
    #                   poped_db)
    # sample_size=1000
    # data = apply(bpop_vals_ed_ln,1,
    #               function(par) rlnorm(sample_size,log(par[2]^2/sqrt(par[3]+par[2]^2)),
    #                                    sqrt(log(par[3]/par[2]^2+1))))
    # colMeans(data)
    # var(data)
    
    
    if bLHS == 0: #Random Sampling
        for k in range(0,sample_size):
            np = size(par)[0]
            if np != 0:
                n = numpy.random.randn(np,1) # normal mean=0 and sd=1
                u = numpy.random.rand(np,1)*2-1 # uniform from -1 to 1
                t = par[:,0] # type of distribution
                c2 = par[:,2] # variance or range of distribution
                
                bUserSpecifiedDistribution = (sum(t==3)>=1) #If at least one user specified distribution
                ret[k,] = (t==0)*par[:,1] + (t==2)*(par[:,1]+u*c2/2) + (t==1)*(par[:,1]+n*c2^(1/2)) + (t==4)*numpy.exp((numpy.log(par[:,1]^2/numpy.sqrt(par[:,2]+par[:,1]^2)))+n*(numpy.sqrt(numpy.log(par[:,2]/par[:,1]^2+1))))
                #(t==4)*par[,2,drop=F]*exp(n*c2^(1/2))
                if sum(t==5)>0: #Truncated normal
                    for i in range(0,size(par)[0]):
                        if t(i)==5:
                            ret[k,i] = getTruncatedNormal(par[i,2],c2[i])
                
                if bUserSpecifiedDistribution:
                    if len(sample_number) == 0:
                        ret[k,] = feval(user_dist_pointer,ret[k,:],t,k,poped_db)
                    else:
                        ret[k,] = feval(user_dist_pointer,ret[k,:],t,sample_number,poped_db)

    elif nvar != 0: #LHS
        ran = numpy.random.rand(sample_size,nvar)
        # method of Stein
        for j in range(0,nvar):
            idx = numpy.random.random_sample(sample_size)
            P=(idx-ran[:,j])/sample_size 
                  # probability of the cdf
            returnArgs_list = [par[j,1],  #point
                                par[j,1]+norm.ppf(P)*numpy.sqrt(par[j,3]), # normal
                                par[j,1]-par[j,3]/2+P*par[j,3], #uniform
                                ret[:,j], #Do nothing
                                numpy.exp((numpy.log(par[j,2]^2/numpy.sqrt(par[j,3]+par[j,2]^2)))+norm.ppf(P)*(numpy.sqrt(numpy.log(par[j,3]/par[j,2]^2+1))))] #log-normal 
                                #par[j,2]*exp(qnorm(P)*sqrt(par[j,3])) #log-normal
            returnArgs = returnArgs_list[par[j,0]+1]
            
        
        if returnArgs is None:
            raise Exception("Unknown distribution for the inverse probability function used in Latin Hypercube Sampling")
        
        ret[:,j] = returnArgs
        
        bUserSpecifiedDistribution = (sum(par[:,0]==3)>=1) #If at least one user specified distribution
        
        if bUserSpecifiedDistribution is True:
            for k in range(0,sample_size):
                if len(sample_number) == 0:
                    ret[k,] = feval(user_dist_pointer,ret[k,:],par[:,0],k,poped_db)
                else:
                    ret[k,] = feval(user_dist_pointer,ret[k,:],par[:,0],sample_number,poped_db)

    return ret



#sign.2 = function(x){
#  s = x/abs(x)*(x!=0)
#return( s) 
#}


