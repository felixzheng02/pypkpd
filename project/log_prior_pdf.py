"""
## Compute the natural log of the PDF for the parameters in an E-family design
## 
## @inheritParams ed_laplace_ofv
## @inheritParams evaluate.fim
## @inheritParams create.poped.database
## @inheritParams Doptim
## @param alpha  A parameter vector.
## @param return_hessian Should the hessian be returned?
## 
# @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
# @example tests/testthat/examples_fcn_doc/examples_log_prior_pdf.R
# @export
## @keywords internal

## Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
import scipy as sp
from project.zeros import zeros
from project.diag_matlab import diag_matlab

def log_prior_pdf(alpha,bpopdescr,ddescr,return_gradient=False,return_hessian=False):
    #returns the logarithm of the probability density for alpha,given the prior
    #and if required the gradient
    #priordescr=Matrix(c(bpopdescr[bpopdescr[,1]!=0,], ddescr[ddescr[,1]!=0,]),nrow=1,byrow=T)
    priordescr = np.concatenate((bpopdescr,ddescr), axis=0)
    priordescr = priordescr[priordescr[:,0]!=0,:]
    if any(priordescr[:,0] > 4 or priordescr[:,0] == 3):
        raise Exception("Specified prior distribution not supported")
    
    mu = priordescr[:,1]	#means
    sigma = priordescr[:,2] #variances
    n = priordescr[:,0] == 1
    u = priordescr[:,0] == 2
    ln = priordescr[:,0] == 4
    
    #normal
    #temp=(-log(sqrt(2*pi*sigma))-(alpha-mu)^2/(2*sigma))*n
    temp_n = np.log(sp.stats.norm.pdf(alpha[n], mean=mu[n], sd=np.sqrt(sigma[n]))) 
    #p=sum(temp[n])
    p = sum(temp_n)
    
    #uniform
    #temp=log((alpha>=mu-sigma/2 & alpha<=mu+sigma/2)*(1/sigma))
    #p=p+sum(temp[u])
    temp_u = np.log((alpha[u]>=mu[u]-sigma[u]/2 & alpha[u]<=mu[u]+sigma[u]/2)*(1/sigma[u]))
    p = p + sum(temp_u)
    
    #log normal
    #temp=log(1/alpha*own_normpdf(-log(mu)+log(alpha),0,sqrt(sigma)))
    #p=p+sum(temp[ln])
    #temp_ln=log(1/alpha[ln]*own_normpdf(-log(mu[ln])+log(alpha[ln]),0,sqrt(sigma[ln])))
    temp_ln = np.log(sp.stats.lognorm.pdf(alpha[ln],np.log(mu[ln]^2/np.sqrt(sigma[ln]+mu[ln]^2)),np.sqrt(np.log(sigma[ln]/mu[ln]^2+1))))
    p = p + sum(temp_ln)
    
    grad = None
    if return_gradient is True: 
        grad = zeros(alpha.size, 1)
        #normal
        grad[n] = -(alpha[n]-mu[n])/sigma[n]
        #uniform
        grad[u] = 0
        #log normal
        temp = -(sigma+np.log(alpha)-np.log(mu))/(alpha*sigma)
        grad[ln] = temp[ln]


    hess = None
    if return_hessian is True:
        diaghess = zeros(alpha.size, 1)
        diaghess[n] = -1/sigma[n]
        diaghess[u] = 0
        temp = (-1+sigma+np.log(alpha)-np.log(mu))/(alpha^2*sigma)
        diaghess[ln] = temp[ln]
        hess = diag_matlab(diaghess)

    
    #return(list( p= p,grad=grad,hess=hess)) 
    ret_args = "p"
    if grad is not None or hess is not None: 
        ret_args = '{"p": p'
    if grad is not None: 
        ret_args = ret_args + ', "grad": grad'
    if hess is not None: 
        ret_args = ret_args + ', "hess": hess'
    if grad is not None or hess is not None: 
        ret_args = ret_args + "}"
            
    return exec("text = " + ret_args)
  


# own_normpdf = function(x,mu,sigma){
#   y = exp(-0.5 * ((x - mu)/sigma)^2) / (sqrt(2*pi) * sigma)
#   return( y ) 
# }
