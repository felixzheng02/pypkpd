"""
create.poped.database.R
Create a PopED database
This function takes the input file (a previously created poped database) supplied by the user, or function arguments, 
and creates a database that can then be used to 
run all other PopED functions.  The function supplies default values to elements of the 
database that are not specified in the
input file or as function arguments. Default arguments are supplied in the Usage section 
(easiest to use a text search to find values you are interested in).  
@inheritParams create_design_space
@param popedInput A PopED database file or an empty list \code{list()}.  List elements should match the values seen in 
the Usage section (the defaults to function arguments). 
@param ff_file  \itemize{
\item \bold{******START OF MODEL DEFINITION OPTIONS**********}
}
A string giving the function name or filename and path of the structural model. 
The filename and the function name must be the same if giving a filename. 
e.g. \code{"ff.PK.1.comp.oral.md.KE"}
@param ff_fun Function describing the structural model. e.g. \code{ff.PK.1.comp.oral.md.KE}. 
@param fg_file A string giving the function name or filename and path of the 
parameter model. 
The filename and the function name must be the same if giving a filename. 
e.g. \code{"parameter.model"}
@param fg_fun Function describing the parameter model. e.g. \code{parameter.model}.
@param fError_file A string giving the function name or filename and path of the 
residual error model. 
The filename and the function name must be the same if giving a filename. 
e.g. \code{"feps.prop"}.
@param fError_fun Function describing the residual error model. e.g. \code{feps.prop}.
#'
@param optsw  \itemize{
\item \bold{******WHAT TO OPTIMIZE**********}}
 Row vector of optimization tasks (1=True,0=False) in the following order: 
(Samples per subject, Sampling schedule, Discrete design variable, Continuous design variable, Number of id per group). 
All elements set to zero => only calculate the FIM with current design
@param xt  \itemize{
\item \bold{******START OF INITIAL DESIGN OPTIONS**********}}
 Matrix defining the initial sampling schedule. 
 Each row is a group/individual.
 If only one vector is supplied, e.g. \code{c(1,2,3,4)}, then all groups will 
have the same initial design. 
@param m Number of groups in the study.  Each individual in a group will have the same design. 
@param x A matrix defining the initial discrete values for the model 
Each row is a group/individual.
@param nx Number of discrete design variables.
@param a Matrix defining the initial continuous covariate values. 
n_rows=number of groups, n_cols=number of covariates.
If the number of rows is one and the number of groups > 1 then all groups are assigned the 
same values.
# @param na The number of covariates in the model.
@param groupsize Vector defining the size of the different groups (num individuals in each group).
If only one number then the number will be the same in every group.
@param ni Vector defining the number of samples for each group. 
@param model_switch Matrix defining which response a certain sampling time belongs to.
@param maxni  \itemize{
\item \bold{******START OF DESIGN SPACE OPTIONS**********}}
Max number of samples per group/individual
@param minni Min number of samples per group/individual 
@param maxgroupsize Vector defining the max size of the different groups (max number of individuals in each group)
@param mingroupsize Vector defining the min size of the different groups (min num individuals in each group) --
@param maxtotgroupsize The total maximal groupsize over all groups
@param mintotgroupsize The total minimal groupsize over all groups
@param maxxt Matrix or single value defining the maximum value for each xt sample.  If a single value is 
supplied then all xt values are given the same maximum value.
@param minxt Matrix or single value defining the minimum value for each xt sample.  If a single value is 
supplied then all xt values are given the same minimum value
@param discrete_x Cell array defining the discrete variables for each x value. 
  See examples in \code{\link{create_design_space}}.
@param discrete_xt Cell array \code{\link{cell}} defining the discrete variables allowed for each xt value.
  Can also be a list of values \code{list(1:10)} (same values allowed for all xt), or a list of lists 
 \code{list(1:10, 2:23, 4:6)} (one for each value in xt). See examples in \code{\link{create_design_space}}.
@param discrete_a Cell array \code{\link{cell}} defining the discrete variables allowed for each a value.
  Can also be a list of values \code{list(1:10)} (same values allowed for all a), or a list of lists 
 \code{list(1:10, 2:23, 4:6)} (one for each value in a). See examples in \code{\link{create_design_space}}.
@param maxa Vector defining the max value for each covariate. If a single value is supplied then
 all a values are given the same max value
@param mina Vector defining the min value for each covariate. If a single value is supplied then
 all a values are given the same max value
@param bUseGrouped_xt Use grouped time points (1=True, 0=False).
@param G_xt Matrix defining the grouping of sample points. Matching integers mean that the points are matched.
@param bUseGrouped_a Use grouped covariates (1=True, 0=False)
@param G_a Matrix defining the grouping of covariates. Matching integers mean that the points are matched.
@param bUseGrouped_x Use grouped discrete design variables (1=True, 0=False).
@param G_x  Matrix defining the grouping of discrete design variables. Matching integers mean that the points are matched.
@param iFIMCalculationType  \itemize{
\item \bold{******START OF FIM CALCULATION OPTIONS**********}}
Fisher Information Matrix type
\itemize{
\item 0=Full FIM
\item 1=Reduced FIM
\item 2=weighted models
\item 3=Loc models
\item 4=reduced FIM with derivative of SD of sigma as in PFIM
\item 5=FULL FIM parameterized with A,B,C matrices & derivative of variance
\item 6=Calculate one model switch at a time, good for large matrices
\item 7=Reduced FIM parameterized with A,B,C matrices & derivative of variance
}
@param iApproximationMethod Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI
@param iFOCENumInd Num individuals in each step of FOCE 
@param prior_fim The prior FIM (added to calculated FIM) 
@param strAutoCorrelationFile Filename and path, or function name, for the Autocorrelation function, 
empty string means no autocorrelation.
@param d_switch  \itemize{
\item \bold{******START OF CRITERION SPECIFICATION OPTIONS**********}}
D-family design (1) or ED-family design (0) (with or without parameter uncertainty) 
@param ofv_calc_type  OFV calculation type for FIM 
\itemize{ 
\item 1 = "D-optimality". Determinant of the FIM: det(FIM)
\item 2 = "A-optimality".  Inverse of the sum of the expected parameter variances: 
1/trace_matrix(inv(FIM)) 
\item 4 = "lnD-optimality".  Natural logarithm of the determinant of the FIM: log(det(FIM)) 
\item 6 = "Ds-optimality". Ratio of the Determinant of the FIM and the Determinant of the uninteresting
rows and columns of the FIM: det(FIM)/det(FIM_u)
\item 7 = Inverse of the sum of the expected parameter RSE: 1/sum(get_rse(FIM,poped_db,use_percent=False))
}
@param ds_index Ds_index is a vector set to 1 if a parameter is uninteresting, otherwise 0.
size=(1,num unfixed parameters). First unfixed bpop, then unfixed d, then unfixed docc and last unfixed sigma. 
Default is the fixed effects being important, everything else not important.  Used in conjunction with
\code{ofv_calc_type=6}.
@param strEDPenaltyFile Penalty function name or path and filename, empty string means no penalty.  
User defined criterion can be defined this way.
@param ofv_fun User defined function used to compute the objective function. The function must have a poped database object as its first
argument and have "..." in its argument list.  Can be referenced as a function or as a file name where the function defined in the file has the same name as the file.
e.g. "cost.txt" has a function named "cost" in it.
@param iEDCalculationType  \itemize{
\item \bold{******START OF E-FAMILY CRITERION SPECIFICATION OPTIONS**********}}
ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
@param ED_samp_size Sample size for E-family sampling 
@param bLHS How to sample from distributions in E-family calculations. 0=Random Sampling, 1=LatinHyperCube --
@param strUserDistributionFile Filename and path, or function name, for user defined distributions for E-family designs 
@param nbpop  \itemize{
\item \bold{******START OF Model parameters  SPECIFICATION OPTIONS**********}}
Number of typical values 
@param NumRanEff Number of IIV parameters. Typically can be computed from other values and not supplied. 
@param NumDocc Number of IOV variance parameters. Typically can be computed from other values and not supplied. 
@param NumOcc Number of occasions. Typically can be computed from other values and not supplied. 
# @param ng The length of the g parameter vector. Typically can be computed from other values and not supplied.
@param bpop Matrix defining the fixed effects, per row (row number = parameter_number) we should have:
\itemize{
\item column 1 the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
 3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal)
\item column 2  defines the mean.
\item column 3 defines the variance of the distribution (or length of uniform distribution).
}
Can also just supply the parameter values as a vector \code{c()} if no uncertainty around the 
parameter value is to be used. The parameter order of  'bpop' is defined in the 'fg_fun' or 'fg_file'. If you use named 
arguments in 'bpop' then the order will be worked out automatically.
@param d Matrix defining the diagonals of the IIV (same logic as for the fixed effects 
matrix bpop to define uncertainty). One can also just supply the parameter values as a \code{c()}. 
The parameter order of 'd' is defined in the 'fg_fun' or 'fg_file'. If you use named 
arguments in 'd' then the order will be worked out automatically.
@param covd Column major vector defining the covariances of the IIV variances. 
That is, from your full IIV matrix  \code{covd =  IIV[lower.tri(IIV)]}. 
@param sigma Matrix defining the variances can covariances of the residual variability terms of the model.
can also just supply the diagonal parameter values (variances) as a \code{c()}. 
@param docc Matrix defining the IOV, the IOV variances and the IOV distribution as for d and bpop. 
@param covdocc Column major vector defining the covariance of the IOV, as in covd. 
@param notfixed_bpop  \itemize{
\item \bold{******START OF Model parameters fixed or not  SPECIFICATION OPTIONS**********}}
Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed). 
The parameter order of 'notfixed_bpop' is defined in the 'fg_fun' or 'fg_file'. If you use named 
arguments in 'notfixed_bpop' then the order will be worked out automatically.
@param notfixed_d Vector defining if a IIV is fixed or not (1=not fixed, 0=fixed). 
The parameter order of 'notfixed_d' is defined in the 'fg_fun' or 'fg_file'. If you use named 
arguments in 'notfixed_d' then the order will be worked out automatically. 
@param notfixed_covd Vector defining if a covariance IIV is fixed or not (1=not fixed, 0=fixed)
@param notfixed_docc Vector defining if an IOV variance is fixed or not (1=not fixed, 0=fixed)  
@param notfixed_covdocc Vector row major order for lower triangular matrix defining if a covariance IOV is fixed or not (1=not fixed, 0=fixed) 
@param notfixed_sigma Vector defining if a residual error parameter is fixed or not (1=not fixed, 0=fixed) 
@param notfixed_covsigma Vector defining if a covariance residual error parameter is fixed or not (1=not fixed, 0=fixed). 
Default is fixed.
@param bUseRandomSearch  \itemize{
\item \bold{******START OF Optimization algorithm  SPECIFICATION OPTIONS**********}}
Use random search (1=True, 0=False)
@param bUseStochasticGradient Use Stochastic Gradient search (1=True, 0=False) 
@param bUseLineSearch Use Line search (1=True, 0=False) 
@param bUseExchangeAlgorithm Use Exchange algorithm (1=True, 0=False)        
@param bUseBFGSMinimizer Use BFGS Minimizer (1=True, 0=False) 
@param EACriteria Exchange Algorithm Criteria, 1 = Modified, 2 = Fedorov 
@param strRunFile Filename and path, or function name, for a run file that is used instead of the regular PopED call. 
@param poped_version  \itemize{
\item \bold{******START OF Labeling and file names  SPECIFICATION OPTIONS**********}}
The current PopED version 
@param modtit The model title 
@param output_file Filename and path of the output file during search 
@param output_function_file Filename suffix of the result function file 
@param strIterationFileName Filename and path for storage of current optimal design 
@param user_data  \itemize{
\item \bold{******START OF Miscellaneous SPECIFICATION OPTIONS**********}}
User defined data structure that, for example could be used to send in data to the model 
@param ourzero Value to interpret as zero in design 
@param dSeed The seed number used for optimization and sampling -- integer or -1 which creates a random seed \code{as.integer(Sys.time())} or None.
@param line_opta Vector for line search on continuous design variables (1=True,0=False)
@param line_optx Vector for line search on discrete design variables (1=True,0=False) 
@param bShowGraphs Use graph output during search
@param use_logfile If a log file should be used (0=False, 1=True)
@param m1_switch Method used to calculate M1 
(0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) 
@param m2_switch Method used to calculate M2
(0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) 
@param hle_switch Method used to calculate linearization of residual error
(0=Complex difference, 1=Central difference, 30=Automatic differentiation) 
@param gradff_switch Method used to calculate the gradient of the model
(0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) 
@param gradfg_switch Method used to calculate the gradient of the parameter vector g
(0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) 
@param grad_all_switch Method used to calculate all the gradients
(0=Complex difference, 1=Central difference) 
@param rsit_output Number of iterations in random search between screen output 
@param sgit_output Number of iterations in stochastic gradient search between screen output 
@param hm1 Step length of derivative of linearized model w.r.t. typical values 
@param hlf Step length of derivative of model w.r.t. g 
@param hlg Step length of derivative of g w.r.t. b
@param hm2 Step length of derivative of variance w.r.t. typical values 
@param hgd Step length of derivative of OFV w.r.t. time 
@param hle Step length of derivative of model w.r.t. sigma
@param AbsTol The absolute tolerance for the diff equation solver
@param RelTol The relative tolerance for the diff equation solver
@param iDiffSolverMethod The diff equation solver method, None as default.
@param bUseMemorySolver If the differential equation results should be stored in memory (1) or not (0)
@param rsit Number of Random search iterations 
@param sgit Number of stochastic gradient iterations
@param intrsit Number of Random search iterations with discrete optimization.
@param intsgit Number of Stochastic Gradient search iterations with discrete optimization 
@param maxrsnullit Iterations until adaptive narrowing in random search
@param convergence_eps Stochastic Gradient convergence value,
(difference in OFV for D-optimal, difference in gradient for ED-optimal)
@param rslxt Random search locality factor for sample times 
@param rsla Random search locality factor for covariates 
@param cfaxt Stochastic Gradient search first step factor for sample times 
@param cfaa Stochastic Gradient search first step factor for covariates 
@param bGreedyGroupOpt Use greedy algorithm for group assignment optimization 
@param EAStepSize Exchange Algorithm StepSize 
@param EANumPoints Exchange Algorithm NumPoints 
@param EAConvergenceCriteria Exchange Algorithm Convergence Limit/Criteria 
@param bEANoReplicates Avoid replicate samples when using Exchange Algorithm 
@param BFGSConvergenceCriteriaMinStep BFGS Minimizer Convergence Criteria Minimum Step 
@param BFGSProjectedGradientTol BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance 
@param BFGSTolerancef BFGS Minimizer Line Search Tolerance f 
@param BFGSToleranceg BFGS Minimizer Line Search Tolerance g 
@param BFGSTolerancex BFGS Minimizer Line Search Tolerance x 
@param ED_diff_it Number of iterations in ED-optimal design to calculate convergence criteria 
@param ED_diff_percent ED-optimal design convergence criteria in percent 
@param line_search_it Number of grid points in the line search 
@param Doptim_iter Number of iterations of full Random search and full Stochastic Gradient if line search is not used 
@param iCompileOption \bold{******START OF PARALLEL OPTIONS**********} Compile options for PopED
\itemize{
\item -1 = No compilation,
\item 0 or 3 = Full compilation,
\item 1 or 4 = Only using MCC (shared lib),
\item 2 or 5 = Only MPI,
\item Option 0,1,2 runs PopED and option 3,4,5 stops after compilation
}
@param iUseParallelMethod Parallel method to use (0 = Matlab PCT, 1 = MPI) 
@param MCC_Dep Additional dependencies used in MCC compilation (mat-files), if several space separated 
@param strExecuteName Compilation output executable name 
@param iNumProcesses Number of processes to use when running in parallel (e.g. 3 = 2 workers, 1 job manager) 
@param iNumChunkDesignEvals Number of design evaluations that should be evaluated in each process before getting new work from job manager
# @param strMatFileInputPrefix The prefix of the input mat file to communicate with the executable 
@param Mat_Out_Pre The prefix of the output mat file to communicate with the executable 
@param strExtraRunOptions Extra options send to e$g. the MPI executable or a batch script, see execute_parallel$m for more information and options 
@param dPollResultTime Polling time to check if the parallel execution is finished 
@param strFunctionInputName The file containing the popedInput structure that should be used to evaluate the designs 
@param bParallelRS If the random search is going to be executed in parallel
@param bParallelSG If the stochastic gradient search is going to be executed in parallel 
@param bParallelMFEA If the modified exchange algorithm is going to be executed in parallel 
@param bParallelLS If the line search is going to be executed in parallel 
@return A PopED database
Author: Caiya Zhang, Yuchen Zheng
"""


# from am.poped_choose import am.poped_choose
import random
import datetime
import numpy as np
import pandas as pd
from os import name
from project.sfg import sfg
from project.cell import cell
from project.ones import ones
from project.size import size
from matpy.num import num
from matpy.matrix import matrix
from project.zeros import zeros
from project.feval import feval
from project.pargen import pargen
from project.util import is_not_none
from project.getfulld import getfulld
from project.fileparts import fileparts
from project.poped_choose import poped_choose
from project.param_choose import param_choose
from project.test_mat_size import test_mat_size
from project.create_design import create_design
from project.create_design_space import create_design_space
from project.find_largest_index import find_largest_index
from project.convert_variables import convert_variables
from project.get_all_params import get_all_params
from project.get_unfixed_params import get_unfixed_params


def reorder_vec(your_vec: matrix, name_order):
    if your_vec.get_datanam() is not None:
        if all(your_vec.get_datanam()) in name_order:
            your_vec = your_vec[name_order[name_order in your_vec.get_datanam()]]

    return your_vec


def create_poped_database(popedInput={}, **kwargs):
    

    
    
    param = list(kwargs.keys())

    
    # for i in range(0, len(param)):
    #     k = param[i]
    #     v = kwargs[k]
    #     exec("%s = %s" % (k, "v"))
    # parameters

    # --------------------------
    # ---- Model definition
    # --------------------------
    keys = list(popedInput.keys())
    # -- Filname and path of the model file --
    ff_file = param_choose(kwargs, None, 0, param, kwargs, "ff_file")

    # if "ff_file" not in param:
    #     ff_file = None
    ff_fun = param_choose(popedInput, None, 0, param, kwargs, "model", "ff_pointer")

    # -- Filname and path of the g parameter file --
    fg_file = param_choose(kwargs, None, 0, param, kwargs, "fg_file")
    # if "fg_file" not in param:
    #     fg_file = None
    fg_fun = param_choose(popedInput, None, 0, param, kwargs, "model", "fg_pointer")
    # -- Filname and path of the error model file --
    fError_file = param_choose(kwargs, None, 0, param, kwargs, "fError_file")
    # if "fError_file" not in param:
    #     fError_file = None
    fError_fun = param_choose(popedInput, None, 0, param, kwargs, "model", "ferror_pointer")

    # --------------------------
    # ---- What to optimize
    # --------------------------

    # -- Vector of optimization tasks (1=True,0=False)
    # (Samples per subject, Sampling schedule, Discrete design variable, Continuous design variable, Number of id per group)
    # -- All elements set to zero => only calculate the FIM with current design --
    optsw = param_choose(popedInput, matrix(np.array([0, 0, 0, 0, 0])), 0, param, kwargs, 'settings', 'optsw')

    # --------------------------
    # ---- Initial Design
    # --------------------------

    # -- Matrix defining the initial sampling schedule --
    
    xt = param_choose(popedInput, "'xt' needs to be defined", 1, param, kwargs, 'design', 'xt')
    # -- Number of groups/individuals --
    # size(,1) find number of rows
    #m=poped_choose(popedInput['m'], size(xt,1)),
    m = param_choose(popedInput, None, 0, param, kwargs, 'design', 'm')
    # -- Matrix defining the initial discrete values --
    #x=poped_choose(popedInput['design']['x'], zeros(m,0)),
    x = param_choose(popedInput, None, 0, param, kwargs, 'design', 'x')
    # -- Number of discrete design variables --
    # size(,1) find number of cols
    #nx=poped_choose(popedInput['nx'], size(x,2)),
    nx = param_choose(popedInput, None, 0, param, kwargs, 'design', 'nx')
    # -- Vector defining the initial covariate values --
    # a=poped_choose(popedInput$design["a"],zeros(m,0)),
    a = param_choose(popedInput, None, 0, param, kwargs, 'design', 'a')
    # number of continuous design variables that are not time (e.g. continuous covariates)
    # na=poped_choose(popedInput['na'],size(a,2)),
    # na=poped_choose(popedInput['na'],None),
    # -- Vector defining the size of the different groups (num individuals in each group) --

    groupsize = param_choose(popedInput, "'groupsize' needs to be defined", 1, param, kwargs, 'design', 'groupsize')
    # -- Vector defining the number of samples for each group --
    # ni=poped_choose(popedInput["design"]ni,matrix(size(xt,2),m,1)),
    ni = param_choose(popedInput, None, 0, param, kwargs, 'design', 'ni')
    # -- Vector defining which response a certain sampling time belongs to --
    # model_switch=poped_choose(popedInput["design"]model_switch,ones(size(xt,1),size(xt,2))),
    model_switch = param_choose(popedInput, None, 0, param, kwargs, 'design', 'model_switch')

    # --------------------------
    # ---- Design space
    # --------------------------

    # -- Max number of samples per group/individual --
    maxni = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'maxni')
    # -- Min number of samples per group/individual --
    minni = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'minni')
    maxtotni = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'maxtotni')
    mintotni = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'mintotni')
    # -- Vector defining the max size of the different groups (max num individuals in each group) --
    maxgroupsize = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'maxgroupsize')
    # -- Vector defining the min size of the different groups (min num individuals in each group) --
    # mingroupsize=poped_choose(popedInput["design"]mingroupsize,ones(m,1)),
    mingroupsize = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'mingroupsize')
    # -- The total maximal groupsize over all groups--
    maxtotgroupsize = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'maxtotgroupsize')
    # -- The total minimal groupsize over all groups--
    mintotgroupsize = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'mintotgroupsize')
    # -- Matrix defining the max value for each sample --
    maxxt = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'maxxt')
    # -- Matrix defining the min value for each sample --
    minxt = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'minxt')
    discrete_xt = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'xt_space')
    # -- Cell defining the discrete variables --
    # discrete_x=poped_choose(popedInput["design"]discrete_x,cell(m,nx)),
    discrete_x = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'discrete_x')
    # -- Vector defining the max value for each covariate --
    maxa = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'maxa')
    # -- Vector defining the min value for each covariate --
    mina = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'mina')
    discrete_a = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'a_space')
    # -- Use grouped time points (1=True, 0=False) --
    bUseGrouped_xt = param_choose(popedInput, False, 0, param, kwargs, 'design_space', 'bUseGrouped_xt')
    # -- Matrix defining the grouping of sample points --
    G_xt = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'G_xt')
    # -- Use grouped covariates (1=True, 0=False) --
    bUseGrouped_a = param_choose(popedInput, False, 0, param, kwargs, 'design_space', 'bUseGrouped_a')
    # -- Matrix defining the grouping of covariates --
    G_a = param_choose(popedInput, None, 0, param, kwargs, 'design_space', 'G_a')
    # -- Use grouped discrete design variables (1=True, 0=False) --
    bUseGrouped_x = param_choose(popedInput, False, 0, param, kwargs, 'design_space', 'bUseGrouped_x')
    # -- Matrix defining the grouping of discrete design variables --
    G_x = param_choose(popedInput, None, 0, param, kwargs, 'design_space', "G_x")

    # --------------------------
    # ---- FIM calculation
    # --------------------------

    # -- Fisher Information Matrix type
    # (0=Full FIM,
    # 1=Reduced FIM,
    # 2=weighted models,
    # 3=Loc models,
    # 4=reduced FIM with derivative of SD of sigma as pfim,
    # 5=FULL FIM parameterized with A,B,C matrices & derivative of variance,
    # 6=Calculate one model switch at a time, good for large matrices,
    # 7=Reduced FIM parameterized with A,B,C matrices & derivative of variance) --
    iFIMCalculationType = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'iFIMCalculationType')
    # -- Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI --
    iApproximationMethod = param_choose(popedInput, 0, 0, param, kwargs, 'settings', 'iApproximationMethod')
    # -- Num individuals in each step of FOCE --
    iFOCENumInd = param_choose(popedInput, 1000, 0, param, kwargs, 'settings', 'iFOCENumInd')
    # -- The prior FIM (added to calculated FIM) --
    prior_fim = param_choose(popedInput, matrix(np.array([0])), 0, param, kwargs, 'settings', 'prior_fim')
    # -- Filname and path for the Autocorrelation function, empty string means no autocorrelation --
    strAutoCorrelationFile = param_choose(popedInput, "", 0, param, kwargs, "model", 'auto_pointer')

    # --------------------------
    # ---- Criterion specification
    # --------------------------

    # -- D-family design (1) or ED-family design (0) (with or without parameter uncertainty) --
    d_switch = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'd_switch')
    # -- OFV calculation type for FIM (1=Determinant of FIM,4=log determinant of FIM,6=determinant of interesting part of FIM (Ds)) --
    ofv_calc_type = param_choose(popedInput, 4, 0, param, kwargs, 'settings', 'ofv_calc_type')
    # -- Ds_index, set index to 1 if a parameter is uninteresting, otherwise 0.
    # size=(1,num unfixed parameters). First unfixed bpop, then unfixed d, then unfixed docc and last unfixed sigma --
    # default is the fixed effects being important
    ds_index = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'ds_index')
    # -- Penalty function, empty string means no penalty.  User defined criterion --
    strEDPenaltyFile = param_choose(popedInput, "", 0, param, kwargs, 'settings', 'strEDPenaltyFile')
    ofv_fun = param_choose(popedInput, None, 0, param, kwargs, 'settings', 'ofv_fun')

    # --------------------------
    # ---- E-family Criterion options
    # --------------------------
    # -- ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
    iEDCalculationType = param_choose(popedInput, 0, 0, param, kwargs, 'settings', 'iEDCalculationType')
    # -- Sample size for E-family sampling --
    ED_samp_size = param_choose(popedInput, 45, 0, param, kwargs, 'settings', 'ED_samp_size')
    # -- How to sample from distributions in E-family calculations. 0=Random Sampling, 1=LatinHyperCube --
    bLHS = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'bLHS')
    # -- Filname and path for user defined distributions for E-family designs --
    strUserDistributionFile = param_choose(popedInput, "", 0, param, kwargs, "model", 'user_distribution_pointer')

    # --------------------------
    # ---- Model parameters
    # --------------------------

    # -- Number of typical values --
    nbpop = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'nbpop')
    # -- Number of IIV parameters --
    NumRanEff = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'NumRanEff')
    # -- Number of IOV variance parameters --
    NumDocc = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'NumDocc')
    # -- Number of occassions --
    NumOcc = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'NumOcc')
    # -- The length of the g parameter vector --
    # ng=popedInput["parameters"]ng,

    # -- Matrix defining the fixed effects, per row (row number = parameter_number),
    # the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
    # 3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal).
    # The second column defines the mean.
    # The third column defines the variance of the distribution.
    # can also just supply the parameter values as a c()

    bpop = param_choose(popedInput, 'bpop must be defined', 1, param, kwargs, 'parameters', 'bpop')
    # -- Matrix defining the diagonals of the IIV (same logic as for the fixed effects) --
    # can also just supply the parameter values as a c()
    # -- vector defining the row major lower triangle of the covariances of the IIV variances --
    d = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'd')
    # set to zero if not defined
    covd = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'covd')
    # -- Matrix defining the variances of the residual variability terms --
    # REQUIRED! No defaults given.
    # can also just supply the diagonal values as a c()
    sigma = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'sigma')
    # -- Matrix defining the IOV, the IOV variances and the IOV distribution --
    docc = param_choose(popedInput, matrix(np.array([np.nan, np.nan, np.nan])), 0, param, kwargs, 'parameters', 'docc')
    # -- Matrix defining the covariance of the IOV --
    if np.array_equal(docc.get_data(), np.array([np.nan, np.nan, np.nan]), equal_nan=True):
        tmp = 0
    else:
        tmp = len(docc.get_data()[1])
    covdocc = param_choose(popedInput, zeros(1, tmp*(tmp-1)/2), 0, param, kwargs, 'parameters', 'covdocc')

    # --------------------------
    # ---- Model parameters fixed or not
    # --------------------------
    # -- Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed) --
    notfixed_bpop = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'notfixed_bpop')
    # -- Vector defining if a IIV is fixed or not (1=not fixed, 0=fixed) --
    notfixed_d = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'notfixed_d')
    # -- Vector defining if a covariance IIV is fixed or not (1=not fixed, 0=fixed) --
    notfixed_covd = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'notfixed_covd')
    # -- Vector defining if an IOV variance is fixed or not (1=not fixed, 0=fixed) --
    notfixed_docc = param_choose(popedInput, None, 0, param, kwargs, 'parameters', 'notfixed_docc')
    # -- Vector row major order for lower triangular matrix defining if a covariance IOV is fixed or not (1=not fixed, 0=fixed) --
    notfixed_covdocc = param_choose(popedInput, zeros(1, len(covdocc)), 0, param, kwargs, 'parameters', 'notfixed_covdocc')
    # -- Vector defining if a residual error parameter is fixed or not (1=not fixed, 0=fixed) --
    notfixed_sigma = param_choose(popedInput, matrix(np.ones(size(sigma)[1])), 0, param, kwargs, 'parameters', 'notfixed_sigma')
    # -- Vector defining if a covariance residual error parameter is fixed or not (1=not fixed, 0=fixed) --
    ## default is fixed
    notfixed_covsigma = param_choose(popedInput, zeros(1, len(notfixed_sigma)*(len(notfixed_sigma)-1)/2), 0, param, kwargs, 'parameters', 'notfixed_covsigma')

    # --------------------------
    # ---- Optimization algorithm choices
    # --------------------------

    # -- Use random search (1=True, 0=False) --
    bUseRandomSearch = param_choose(popedInput, True, 0, param, kwargs, 'settings', 'bUseRandomSearch')
    # -- Use Stochastic Gradient search (1=True, 0=False) --
    bUseStochasticGradient = param_choose(popedInput, True, 0, param, kwargs, 'settings', 'bUseStochasticGradient')
    # -- Use Line search (1=True, 0=False) --
    bUseLineSearch = param_choose(popedInput, True, 0, param, kwargs, 'settings', 'bUseLineSearch')
    # -- Use Exchange algorithm (1=True, 0=False) --
    bUseExchangeAlgorithm = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'bUseExchangeAlgorithm')
    # -- Use BFGS Minimizer (1=True, 0=False) --
    bUseBFGSMinimizer = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'bUseBFGSMinimizer')
    # -- Exchange Algorithm Criteria, 1 = Modified, 2 = Fedorov --
    EACriteria = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'EACriteria')
    # -- Filename and path for a run file that is used instead of the regular PopED call --
    strRunFile = param_choose(popedInput, "", 0, param, kwargs, 'settings', 'run_file_pointer')

    # --------------------------
    # ---- Labeling and file names
    # --------------------------

    # -- The current PopED version --
# ！！！没写 packageVersion("PopED")
    poped_version = param_choose(popedInput, "0.0.2", 0, param, kwargs, 'settings', 'poped_version')
    # -- The model title --
    modtit = param_choose(popedInput, 'PopED model', 0, param, kwargs, 'settings', 'modtit')
    # -- Filname and path of the output file during search --
    output_file = param_choose(popedInput, "PopED_output_summary", 0, param, kwargs, 'settings', 'output_file')
    # -- Filname suffix of the result function file --
    output_function_file = param_choose(popedInput, "PopED_output_", 0, param, kwargs, 'settings', 'output_function_file')
    # -- Filename and path for storage of current optimal design --
    strIterationFileName = param_choose(popedInput, "PopED_current.R", 0, param, kwargs, 'settings', 'strIterationFileName')

    # --------------------------
    # ---- Misc options
    # --------------------------
    # -- User defined data structure that, for example could be used to send in data to the model --
    user_data = param_choose(popedInput, cell(0, 0), 0, param, kwargs, 'settings', 'user_data')
    # -- Value to interpret as zero in design --
    ourzero = param_choose(popedInput, 1e-5, 0, param, kwargs, 'settings', 'ourzero')
    # ourzero=param_choose(popedInput$ourzero,0),
    # -- The seed number used for optimization and sampling -- integer or -1 which creates a random seed
    dSeed = param_choose(popedInput, None, 0, param, kwargs, 'settings', 'dSeed')
    # -- Vector for line search on continuous design variables (1=True,0=False) --
    line_opta = param_choose(popedInput, None, 0, param, kwargs, 'settings', 'line_opta')
    # -- Vector for line search on discrete design variables (1=True,0=False) --
    line_optx = param_choose(popedInput, None, 0, param, kwargs, 'settings', 'line_optx')  # matrix(0,0,1)
    # -- Use graph output during search --
    bShowGraphs = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'bShowGraphs')
    # -- If a log file should be used (0=False, 1=True) --
    use_logfile = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'use_logfile')
    # -- Method used to calculate M1
    # (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
    m1_switch = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'm1_switch')
    # -- Method used to calculate M2
    # (0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
    m2_switch = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'm2_switch')
    # -- Method used to calculate linearization of residual error
    # (0=Complex difference, 1=Central difference, 30=Automatic differentiation) --
    hle_switch = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'hle_switch')
    # -- Method used to calculate the gradient of the model
    # (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
    gradff_switch = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'gradff_switch')
    # -- Method used to calculate the gradient of the parameter vector g
    # (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
    gradfg_switch = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'gradfg_switch')
    # -- Method used to calculate all the gradients
    # (0=Complex difference, 1=Central difference) --
    grad_all_switch = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'grad_all_switch')
    # -- Number of iterations in random search between screen output --
    rsit_output = param_choose(popedInput, 5, 0, param, kwargs, 'settings', 'rsit_output')
    # -- Number of iterations in stochastic gradient search between screen output --
    sgit_output = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'sgit_output')
    # -- Step length of derivative of linearized model w.r.t. typical values --
    hm1 = param_choose(popedInput, 0.00001, 0, param, kwargs, 'settings', 'hm1')
    # -- Step length of derivative of model w.r.t. g --
    hlf = param_choose(popedInput, 0.00001, 0, param, kwargs, 'settings', 'hlf')
    # -- Step length of derivative of g w.r.t. b --
    hlg = param_choose(popedInput, 0.00001, 0, param, kwargs, 'settings', 'hlg')
    # -- Step length of derivative of variance w.r.t. typical values --
    hm2 = param_choose(popedInput, 0.00001, 0, param, kwargs, 'settings', 'hm2')
    # -- Step length of derivative of OFV w.r.t. time --
    hgd = param_choose(popedInput, 0.00001, 0, param, kwargs, 'settings', 'hgd')
    # -- Step length of derivative of model w.r.t. sigma --
    hle = param_choose(popedInput, 0.00001, 0, param, kwargs, 'settings', 'hle')
    # -- The absolute tolerance for the diff equation solver --
    AbsTol = param_choose(popedInput, 0.000001, 0, param, kwargs, 'settings', 'AbsTol')
    # -- The relative tolerance for the diff equation solver --
    RelTol = param_choose(popedInput, 0.000001, 0, param, kwargs, 'settings', 'RelTol')
    # -- The diff equation solver method, 0, no other option --
    iDiffSolverMethod = param_choose(popedInput, None, 0, param, kwargs, 'settings', 'iDiffSolverMethod')
    # -- If the differential equation results should be stored in memory (1) or not (0) --
    bUseMemorySolver = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'bUseMemorySolver')
    # -- Number of Random search iterations --
    rsit = param_choose(popedInput, 300, 0, param, kwargs, 'settings', 'rsit')
    # -- Number of Stochastic gradient search iterations --
    sgit = param_choose(popedInput, 150, 0, param, kwargs, 'settings', 'sgit')
    # -- Number of Random search iterations with discrete optimization --
    intrsit = param_choose(popedInput, 250, 0, param, kwargs, 'settings', 'intrsit')
    # -- Number of Stochastic Gradient search iterations with discrete optimization --
    intsgit = param_choose(popedInput, 50, 0, param, kwargs, 'settings', 'intsgit')
    # -- Iterations until adaptive narrowing in random search --
    maxrsnullit = param_choose(popedInput, 50, 0, param, kwargs, 'settings', 'maxrsnullit')
    # -- Stoachstic Gradient convergence value,
    # (difference in OFV for D-optimal, difference in gradient for ED-optimal) --
    convergence_eps = param_choose(popedInput, 1e-08, 0, param, kwargs, 'settings', 'convergence_eps')
    # -- Random search locality factor for sample times --
    rslxt = param_choose(popedInput, 10, 0, param, kwargs, 'settings', 'rslxt')
    # -- Random search locality factor for covariates --
    rsla = param_choose(popedInput, 10, 0, param, kwargs, 'settings', 'rsla')
    # -- Stochastic Gradient search first step factor for sample times --
    cfaxt = param_choose(popedInput, 0.001, 0, param, kwargs, 'settings', 'cfaxt')
    # -- Stochastic Gradient search first step factor for covariates --
    cfaa = param_choose(popedInput, 0.001, 0, param, kwargs, 'settings', 'cfaa')
    # -- Use greedy algorithm for group assignment optimization --
    bGreedyGroupOpt = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'bGreedyGroupOpt')
    # -- Exchange Algorithm StepSize --
    EAStepSize = param_choose(popedInput, 0.01, 0, param, kwargs, 'settings', 'EAStepSize')
    # -- Exchange Algorithm NumPoints --
    EANumPoints = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'EANumPoints')
    # -- Exchange Algorithm Convergence Limit/Criteria --
    EAConvergenceCriteria = param_choose(popedInput, 1e-20, 0, param, kwargs, 'settings', 'EAConvergenceCriteria')
    # -- Avoid replicate samples when using Exchange Algorithm --
    bEANoReplicates = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'bEANoReplicates')
    # -- BFGS Minimizer Convergence Criteria Minimum Step --
    BFGSConvergenceCriteriaMinStep = None,
    # param_choose(popedInput["settings"]BFGSConvergenceCriteriaMinStep,1e-08),
    # -- BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance --
    BFGSProjectedGradientTol = param_choose(popedInput, 0.0001, 0, param, kwargs, 'settings', 'BFGSProjectedGradientTol')
    # -- BFGS Minimizer Line Search Tolerance f --
    BFGSTolerancef = param_choose(popedInput, 0.001, 0, param, kwargs, 'settings', 'BFGSTolerancef')
    # -- BFGS Minimizer Line Search Tolerance g --
    BFGSToleranceg = param_choose(popedInput, 0.9, 0, param, kwargs, 'settings', 'BFGSToleranceg')
    # -- BFGS Minimizer Line Search Tolerance x --
    BFGSTolerancex = param_choose(popedInput, 0.1, 0, param, kwargs, 'settings', 'BFGSTolerancex')
    # -- Number of iterations in ED-optimal design to calculate convergence criteria --
    ED_diff_it = param_choose(popedInput, 30, 0, param, kwargs, 'settings', 'ED_diff_it')
    # -- ED-optimal design convergence criteria in percent --
    ED_diff_percent = param_choose(popedInput, 10, 0, param, kwargs, 'settings', 'ED_diff_percent')
    # -- Number of grid points in the line search --
    line_search_it = param_choose(popedInput, 50, 0, param, kwargs, 'settings', 'ls_step_size')
    # -- Number of iterations of full Random search and full Stochastic Gradient if line search is not used --
    Doptim_iter = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'iNumSearchIterationsIfNotLineSearch')

    # --------------------------
    # -- Parallel options for PopED -- --
    # --------------------------
    #     ## -- Compile option for PopED
    #     ## -1 = No compilation,
    #     ## 0 or 3 = Full compilation,
    #     ## 1 or 4 = Only using MCC (shared lib),
    #     ## 2 or 5 = Only MPI,
    #     ## Option 0,1,2 runs PopED and option 3,4,5 stops after compilation --
    iCompileOption = param_choose(popedInput, -1, 0, param, kwargs, 'settings', 'parallel', 'iCompileOption')
    # -- Parallel method to use (0 = Matlab PCT, 1 = MPI) --
    iUseParallelMethod = param_choose(popedInput, 1, 0, param, kwargs, 'settings', 'parallel', 'iUseParallelMethod')
    # -- Additional dependencies used in MCC compilation (mat-files), if several space separated --
    MCC_Dep = None,
    # param_choose(popedInput["settings"]parallel$strAdditionalMCCCompilerDependencies, ''),
    # -- Compilation output executable name --
    strExecuteName = param_choose(popedInput, 'calc_fim.exe', 0, param, kwargs, 'settings', 'parallel', 'strExecuteName')
    # -- Number of processes to use when running in parallel (e.g. 3 = 2 workers, 1 job manager) --
    iNumProcesses = param_choose(popedInput, 2, 0, param, kwargs, 'settings', 'parallel', 'iNumProcesses')
    # -- Number of design evaluations that should be evaluated in each process before getting new work from job manager --
    iNumChunkDesignEvals = param_choose(popedInput, -2, 0, param, kwargs, 'settings', 'parallel', 'iNumChunkDesignEvals')
    # -- The prefix of the input mat file to communicate with the executable --
    # strMatFileInputPrefix = param_choose(
    #   popedInput["settings"]parallel$strMatFileInputPrefix,
    #   'parallel_input'),
    # -- The prefix of the output mat file to communicate with the executable --
    Mat_Out_Pre = param_choose(popedInput, 'parallel_output', 0, param, kwargs, 'settings', 'parallel', 'strMatFileOutputPrefix')
    # -- Extra options send to e$g. the MPI executable or a batch script, see execute_parallel$m for more information and options --
    strExtraRunOptions = param_choose(popedInput, '', 0, param, kwargs, 'settings', 'parallel', 'strExtraRunOptions')
    # -- Polling time to check if the parallel execution is finished --
    dPollResultTime = param_choose(popedInput, 0.1, 0, param, kwargs, 'settings', 'parallel', 'dPollResultTime')
    # -- The file containing the popedInput structure that should be used to evaluate the designs --
    strFunctionInputName = param_choose(popedInput, 'function_input', 0, param, kwargs, 'settings', 'parallel', 'strFunctionInputName')
    # -- If the random search is going to be executed in parallel --
    bParallelRS = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'parallel', 'bParallelRS')
    # -- If the stochastic gradient search is going to be executed in parallel --
    bParallelSG = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'parallel', 'bParallelSG')
    # -- If the modified exchange algorithm is going to be executed in parallel --
    bParallelMFEA = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'parallel', 'bParallelMFEA')
    # -- If the line search is going to be executed in parallel --
    bParallelLS = param_choose(popedInput, False, 0, param, kwargs, 'settings', 'parallel', 'bParallelLS')

#####-------------- main part --------------#####
    poped_db = {}
# five main headings for database
#     poped_db = list(design=None,
#                      design_space=None,
#                      models=None,
#                      parameters=None,
#                      settings=None)

#     # update popedInput with options supplied in function
#     called_args = match.call()
#     default_args = formals()
#     for(i in names(called_args)[-1]){
#       if(length(grep("^popedInput$",capture.output(default_args[i])))==1) {
#         eval(parse(text=paste(capture.output(default_args[i]),"=",called_args[i])))
#       }
#     }

# modifyList(settings, list()$settings)

# compare to a default input function.
#   ## -- Filname and path of the model file --
#   popedInput$ff_file='ff'
#   ## -- Filname and path of the parameter file --
#   popedInput$fg_file='sfg'
#   ## -- Filname and path of the residual error model file --
#   popedInput$fError_file='feps.add.prop'
#   ## -- The model title --
#   popedInput$modtit='Sigmoidal Emax model'

    poped_db["settings"] = {}
    poped_db["settings"]["poped_version"] = poped_version

    if BFGSConvergenceCriteriaMinStep is None:
        BFGSConvergenceCriteriaMinStep = param_choose(popedInput, 1e-08, 0, "settings", "BFGSConvergenceCriteriaMinStep")

    if MCC_Dep is None:
        MCC_Dep = param_choose(popedInput, "", 0, "settings", "parallel", "strAdditionalMCCCompilerDependencies")

    poped_db["model"] = {}
    poped_db["model"]["user_distribution_pointer"] = ""

    if str(strUserDistributionFile) != "" and strUserDistributionFile is not None:
        if strUserDistributionFile is not None:
            poped_db["model"]["user_distribution_pointer"] = strUserDistributionFile
        else:
            exec(open(strUserDistributionFile).read())
            returnArgs = fileparts(strUserDistributionFile)
            strUserDistFilePath = list(returnArgs.values())[0]
            strUserDistFilename = list(returnArgs.values())[1]
            poped_db["model"]["user_distribution_pointer"] = strUserDistFilename

    # should be removed
    if (np.array(size(x)) == 0).any():
        x = None
        G_x = None
        discrete_x = None

    # should be removed
    if (np.array(size(a)) == 0).any():
        a = None
        G_a = None
        mina = None
        maxa = None

    design = create_design(xt=xt,
                           groupsize=groupsize,
                           m=m,
                           x=x,
                           a=a,
                           ni=ni,
                           model_switch=model_switch)

    design_space = create_design_space(design,
                                       maxni=maxni,
                                       minni=minni,
                                       maxtotni=maxtotni,
                                       mintotni=mintotni,
                                       maxgroupsize=maxgroupsize,
                                       mingroupsize=mingroupsize,
                                       maxtotgroupsize=maxtotgroupsize,
                                       mintotgroupsize=mintotgroupsize,
                                       maxxt=maxxt,
                                       minxt=minxt,
                                       xt_space=discrete_xt,
                                       maxa=maxa,
                                       mina=mina,
                                       a_space=discrete_a,
                                       x_space=discrete_x,
                                       use_grouped_xt=bUseGrouped_xt,
                                       grouped_xt=G_xt,
                                       use_grouped_a=bUseGrouped_a,
                                       grouped_a=G_a,
                                       use_grouped_x=bUseGrouped_x,
                                       grouped_x=G_x,
                                       our_zero=ourzero)

    design = design_space["design"]
    design_space = design_space["design_space"]

    # all of this should be replaced with using the names used in create_design_space as function arguments
    if is_not_none(design_space, "use_grouped_a"):
        design_space["bUseGrouped_a"] = design_space["use_grouped_a"]
    design_space["use_grouped_a"] = None

    if is_not_none(design_space, "use_grouped_x"):
        design_space["bUseGrouped_x"] = design_space["use_grouped_x"]
    design_space["use_grouped_x"] = None

    if is_not_none(design_space, "use_grouped_xt"):
        design_space["bUseGrouped_xt"] = design_space["use_grouped_xt"]
    design_space["use_grouped_xt"] = None

    if is_not_none(design_space, "grouped_a"):
        design_space["G_a"] = design_space["grouped_a"]
    design_space["grouped_a"] = None

    if is_not_none(design_space, "grouped_x"):
        design_space["G_x"] = design_space["grouped_x"]
    design_space["grouped_x"] = None

    if is_not_none(design_space, "grouped_xt"):
        design_space["G_xt"] = design_space["grouped_xt"]
    design_space["grouped_xt"] = None

    if is_not_none(design_space, "x_space"):
        design_space["discrete_x"] = design_space["x_space"]
    design_space["x_space"] = None

# design_space$maxni = max(design_space$maxni)
# design_space$minni = min(design_space$minni)

    # should be removed
    if ~is_not_none(design, "x"):
        design["x"] = zeros(design["m"], 0)
        design_space["G_x"] = design["x"]
        design_space["bUseGrouped_x"] = False
        design_space["discrete_x"] = cell(design["m"], 0)

    # should be removed
    if ~is_not_none(design, "a"):
        design["a"] = zeros(design["m"], 0)
        design_space["G_a"] = design["a"]
        design_space["bUseGrouped_a"] = False
        design_space["mina"] = design["a"]
        design_space["maxa"] = design["a"]

    poped_db["design"] = design
    poped_db["design_space"] = design_space

    # poped_db$m = poped_db["design"]m  # should be removed only in design
# poped_db$nx = param_choose(nx,size(design$x,2)) # should be removed, not needed or in design
# poped_db$na = param_choose(na,size(design$a,2)) # should be removed, not needed or in design
    poped_db["settings"] = {}
    poped_db["settings"]["bLHS"] = bLHS

# poped_db$discrete_x = design_space$discrete_x # should be removed only in design_space

# poped_db$maxni=max(design_space$maxni) # should be only in design_space and called maxmaxni if needed
# poped_db$minni=min(design_space$minni) # should be only in design_space and called minminni if needed

# poped_db$bUseGrouped_xt = design_space$bUseGrouped_xt # should be only in design_space
# poped_db$bUseGrouped_a  = design_space$bUseGrouped_a # should be only in design_space
# poped_db$bUseGrouped_x  = design_space$bUseGrouped_x # should be only in design_space

    poped_db["settings"]["d_switch"] = d_switch
    poped_db["settings"]["iApproximationMethod"] = iApproximationMethod

    poped_db["settings"]["iFOCENumInd"] = iFOCENumInd
    poped_db["settings"]["bUseRandomSearch"] = bUseRandomSearch
    poped_db["settings"]["bUseStochasticGradient"] = bUseStochasticGradient
    poped_db["settings"]["bUseLineSearch"] = bUseLineSearch
    poped_db["settings"]["bUseExchangeAlgorithm"] = bUseExchangeAlgorithm
    poped_db["settings"]["bUseBFGSMinimizer"] = bUseBFGSMinimizer

    poped_db["settings"]["iEDCalculationType"] = iEDCalculationType

    poped_db["settings"]["BFGSConvergenceCriteriaMinStep"] = BFGSConvergenceCriteriaMinStep
    poped_db["settings"]["BFGSProjectedGradientTol"] = BFGSProjectedGradientTol
    poped_db["settings"]["BFGSTolerancef"] = BFGSTolerancef
    poped_db["settings"]["BFGSToleranceg"] = BFGSToleranceg
    poped_db["settings"]["BFGSTolerancex"] = BFGSTolerancex

    poped_db["parameters"] = {}
    poped_db["parameters"]["covdocc"] = covdocc
    poped_db["parameters"]["notfixed_covdocc"] = notfixed_covdocc
    poped_db["parameters"]["notfixed_covsigma"] = notfixed_covsigma

    poped_db["settings"]["parallel"] = {}
    poped_db["settings"]["parallel"]["iCompileOption"] = iCompileOption
    poped_db["settings"]["parallel"]["strAdditionalMCCCompilerDependencies"] = MCC_Dep
    poped_db["settings"]["parallel"]["iUseParallelMethod"] = iUseParallelMethod
    poped_db["settings"]["parallel"]["strExecuteName"] = strExecuteName
    poped_db["settings"]["parallel"]["iNumProcesses"] = iNumProcesses
    poped_db["settings"]["parallel"]["iNumChunkDesignEvals"] = iNumChunkDesignEvals
    #poped_db["settings"]["parallel"]["strMatFileInputPrefix"] = strMatFileInputPrefix
    poped_db["settings"]["parallel"]["strMatFileOutputPrefix"] = Mat_Out_Pre
    poped_db["settings"]["parallel"]["strExtraRunOptions"] = strExtraRunOptions
    poped_db["settings"]["parallel"]["dPollResultTime"] = dPollResultTime
    poped_db["settings"]["parallel"]["strFunctionInputName"] = strFunctionInputName
    poped_db["settings"]["parallel"]["bParallelRS"] = bParallelRS
    poped_db["settings"]["parallel"]["bParallelSG"] = bParallelSG
    poped_db["settings"]["parallel"]["bParallelLS"] = bParallelLS
    poped_db["settings"]["parallel"]["bParallelMFEA"] = bParallelMFEA

    poped_db["settings"]["hm1"] = hm1
    poped_db["settings"]["hlf"] = hlf
    poped_db["settings"]["hlg"] = hlg
    poped_db["settings"]["hm2"] = hm2
    poped_db["settings"]["hgd"] = hgd
    poped_db["settings"]["hle"] = hle

    poped_db["settings"]["AbsTol"] = AbsTol
    poped_db["settings"]["RelTol"] = RelTol
    poped_db["settings"]["iDiffSolverMethod"] = iDiffSolverMethod
    # Temp thing for memory solvers
    poped_db["settings"]["bUseMemorySolver"] = bUseMemorySolver
    poped_db["settings"]["solved_solutions"] = cell(0, 0)
    poped_db["settings"]["maxtime"] = max(max(poped_db["design_space"]["maxxt"])) + str(poped_db["settings"]["hgd"])

    poped_db["settings"]["iFIMCalculationType"] = iFIMCalculationType
    poped_db["settings"]["rsit"] = rsit
    poped_db["settings"]["sgit"] = sgit
    poped_db["settings"]["intrsit"] = intrsit
    poped_db["settings"]["intsgit"] = intsgit
    poped_db["settings"]["maxrsnullit"] = maxrsnullit
    poped_db["settings"]["convergence_eps"] = convergence_eps
    poped_db["settings"]["rslxt"] = rslxt
    poped_db["settings"]["rsla"] = rsla
    poped_db["settings"]["cfaxt"] = cfaxt
    poped_db["settings"]["cfaa"] = cfaa

    poped_db["settings"]["EACriteria"] = EACriteria
    poped_db["settings"]["EAStepSize"] = EAStepSize
    poped_db["settings"]["EANumPoints"] = EANumPoints
    poped_db["settings"]["EAConvergenceCriteria"] = EAConvergenceCriteria

    poped_db["settings"]["ED_samp_size"] = ED_samp_size
    poped_db["settings"]["ED_diff_it"] = ED_diff_it
    poped_db["settings"]["ED_diff_percent"] = ED_diff_percent
    poped_db["settings"]["ls_step_size"] = line_search_it

    poped_db["settings"]["ofv_calc_type"] = ofv_calc_type

    poped_db["settings"]["iNumSearchIterationsIfNotLineSearch"] = Doptim_iter

    poped_db["settings"]["ourzero"] = ourzero
    poped_db["settings"]["rsit_output"] = rsit_output
    poped_db["settings"]["sgit_output"] = sgit_output

    if callable(fg_fun):
        poped_db["model"]["fg_pointer"] = fg_fun
    elif type(fg_fun) is str:
        if fg_fun is not None:
            poped_db["model"]["fg_pointer"] = fg_fun
    else:
        try:
            fg_file
        except NameError:
            # exec(open(fg_file).read())
            returnArgs = fileparts(fg_file)
            strfgModelFilePath = list(returnArgs.values())[0]
            strfgModelFilename = list(returnArgs.values())[1]
            # if (~strcmp(strfgModelFilePath,''))
            # cd(strfgModelFilePath);
            # end
            poped_db["model"]["fg_pointer"] = strfgModelFilename
        else:
            poped_db["model"]["fg_pointer"] = fg_file

    poped_db["settings"]["ed_penalty_pointer"] = zeros(1, 0)
    if str(strEDPenaltyFile) != "":
        if strEDPenaltyFile is not None:
            poped_db["settings"]["ed_penalty_pointer"] = strEDPenaltyFile
        else:
            exec(open(popedInput["strEDPenaltyFile"]).read())
            returnArgs = fileparts(popedInput["strEDPenaltyFile"])
            strEDPenaltyFilePath = list(returnArgs.values())[0]
            strEDPenaltyFilename = list(returnArgs.values())[1]
    # if (~strcmp(strEDPenaltyFilePath,''))
    # cd(strEDPenaltyFilePath);
    # end
            poped_db["settings"]["ed_penalty_pointer"] = strEDPenaltyFilename

    # if(is.null(ofv_fun) or is.function(ofv_fun)){
#   poped_db["settings"]ofv_fun = ofv_fun
# } else {
#   stop("ofv_fun must be a function or None")
# }

    if ofv_fun is None or callable(ofv_fun):
        poped_db["settings"]["ofv_fun"] = ofv_fun
    else:
        # source explicit file
        # here I assume that function in file has same name as filename minus .txt and pathnames
        if ofv_fun is not None:
            exec(open(str(ofv_fun)).read())
            poped_db["settings"]["ofv_fun"] = eval('text=fileparts(ofv_fun)["filename"]')
        else:
            raise Exception("ofv_fun is not a function or None, and no file with that name was found")

    # if(is.function(ofv_fun)){
#   poped_db["settings"]ofv_fun = ofv_fun
# } else if(exists(ofv_fun)){
#   poped_db["settings"]ofv_fun = ofv_fun
# } else {
#   source(ofv_fun)
#   returnArgs =  fileparts(ofv_fun)
#   strffModelFilePath = returnArgs[1]
#   strffModelFilename  = returnArgs[1]
#   ## if (~strcmp(strffModelFilePath,''))
#   ##    cd(strffModelFilePath);
#   ## end
#   poped_db["settings"]ofv_fun = strffModelFilename
# }

    poped_db["model"]["auto_pointer"] = zeros(1, 0)
    #poped_db["model"]["auto_pointer"] = ''
    if strAutoCorrelationFile != "":
        if str(strAutoCorrelationFile) != "":
            if strAutoCorrelationFile is not None:
                poped_db["model"]["auto_pointer"] = strAutoCorrelationFile
            else:
                exec(open(popedInput["strAutoCorrelationFile"]).read())
                returnArgs = fileparts(popedInput["strAutoCorrelationFile"])
                strAutoCorrelationFilePath = list(returnArgs.values())[0]
                strAutoCorrelationFilename = list(returnArgs.values())[1]
            # if (~strcmp(strAutoCorrelationFilePath,''))
            # cd(strAutoCorrelationFilePath);
            # end
                poped_db["model"]["auto_pointer"] = strAutoCorrelationFilename

    if callable(ff_fun):
        poped_db["model"]["ff_pointer"] = ff_fun
    elif type(ff_fun) is str:
        if ff_fun is not None:
            poped_db["model"]["ff_pointer"] = ff_fun
    else:
        try:
            ff_file
        except NameError:
            # exec(open(ff_file).read())
            returnArgs = fileparts(ff_file)
            strffModelFilePath = list(returnArgs.values())[0]
            strffModelFilename = list(returnArgs.values())[1]
            # if (~strcmp(strffModelFilePath,''))
            # cd(strffModelFilePath);
            # end
            poped_db["model"]["ff_pointer"] = strffModelFilename
        else:
            poped_db["model"]["ff_pointer"] = ff_file

    # Check if there is any sub models defined
    for key in popedInput:
        if "SubModels" in popedInput.keys():
            i = 1
            for elem_key in popedInput["SubModels"]:
                if elem_key == ("ff_file%d" % i):
                    # ok<np.nanSGU>
                    exec(eval('popedInput["SubModels"]["ff_file"]%d' % i))
                    returnArgs = fileparts(eval('popedInput["SubModels"]["ff_file"]%d' % i))  # ok<np.nanSGU>
                    strffModelFilePath = list(returnArgs.values())[0]
                    strffModelFilename = list(returnArgs.values())[1]
                    # if (~strcmp(strffModelFilePath,''))
                    # cd(strffModelFilePath);
                    # end
                    poped_db["model"]["subffPointers"]["".join(["ff_pointer", str(i)])] = strffModelFilename
                    i = i + 1

    if callable(fError_fun):
        poped_db["model"]["ferror_pointer"] = fError_fun
    elif type(fError_fun) is str:
        if fError_fun is not None:
            poped_db["model"]["ferror_pointer"] = fError_fun
    else:
        try:
            fError_file
        except NameError:
            # exec(open(fError_file).read())
            returnArgs = fileparts(fError_file)
            strErrorModelFilePath = list(returnArgs.values())[0]
            strErrorModelFilename = list(returnArgs.values())[1]
            # if (~strcmp(strErrorModelFilePath,''))
            # cd(strErrorModelFilePath);
            # end
            poped_db["model"]["ferror_pointer"] = strErrorModelFilename
        else:
            poped_db["model"]["ferror_pointer"] = fError_file

    # %Set the model file string path
# poped_db$model_file = ff_file
##   model_file = eval('functions(poped_db.ff_pointer)');
# if (~strcmp(model_file.file,''))
##       poped_db.model_file = eval('char(model_file.file)');
# else
##       poped_db.model_file = eval('char(model_file.function)');
# end

# ==================================
# Initialize the randomization
# ==================================
    if dSeed is not None:
        if dSeed == -1:
            poped_db["settings"]["dSeed"] = datetime.datetime.now()
        else:
            poped_db["settings"]["dSeed"] = dSeed
        random.seed(poped_db["settings"]["dSeed"])

    poped_db["parameters"]["nbpop"] = poped_choose(nbpop, find_largest_index(poped_db["model"]["fg_pointer"], "bpop"), 0)
    poped_db["parameters"]["NumRanEff"] = poped_choose(NumRanEff, find_largest_index(poped_db["model"]["fg_pointer"], "b"), 0)
    poped_db["parameters"]["NumDocc"] = poped_choose(NumDocc, find_largest_index(poped_db["model"]["fg_pointer"], "bocc", mat=True, mat_row=True), 0)
    poped_db["parameters"]["NumOcc"] = poped_choose(NumOcc, find_largest_index(poped_db["model"]["fg_pointer"], "bocc", mat=True, mat_row=False), 0)
    
    poped_db["parameters"]["ng"] = eval(poped_db["model"]["fg_pointer"] + "(num(0), num(0), num(0), num(0)," + str(zeros(poped_db["parameters"]["NumDocc"], poped_db["parameters"]["NumOcc"])) + ")").get_size()
    
    docc_arr = np.array([1])
    d_arr = np.array([1])
    bpop_arr = np.array([1])
    poped_db["parameters"]["notfixed_docc"] = poped_choose(notfixed_docc, 
        matrix(np.ones([1, param_choose(poped_db, 0, 0, None, None, "parameters", "NumDocc")])), 0)
    poped_db["parameters"]["notfixed_d"] = poped_choose(notfixed_d, 
        matrix(np.ones([1, param_choose(poped_db, 0, 0, None, None, "parameters", "NumRanEff")])), 0)
    poped_db["parameters"]["notfixed_bpop"] = poped_choose(notfixed_bpop, 
        matrix(np.ones([1, param_choose(poped_db, 0, 0, None, None, "parameters", "nbpop")])), 0)


# reorder named values
    fg_names = eval(poped_db["model"]["fg_pointer"] + "(num(1), num(1), num(1), num(1)," + str(ones(param_choose(poped_db, 0, 0, None, None, "parameters", "NumDocc"), param_choose(poped_db, 0, 0, None, None, "parameters", "NumDocc"))) + ")").get_datanam()
    
    poped_db["parameters"]["notfixed_bpop"] = reorder_vec(poped_db["parameters"]["notfixed_bpop"], fg_names)

    poped_db["parameters"]["notfixed_d"] = reorder_vec(poped_db["parameters"]["notfixed_d"], fg_names)

    # we have just the parameter values not the uncertainty
    if size(d)[0] == 1 and size(d)[1] == poped_db["parameters"]["NumRanEff"]:
        d_descr = zeros(poped_db["parameters"]["NumRanEff"], 3)

        # reorder named values
        d = reorder_vec(d, fg_names)
        
        #set_data!!!!!!!!
        d_descr.set_data()[:,1] = d
        d_descr.set_data()[:,0] = 0  # point values
        d_descr.set_data()[:,2] = 0  # variance
        d_descr.set_colnam(d.get_datanam())
        # d_descr.columns.values = d.keys()
        d = d_descr

    # we have just the parameter values not the uncertainty
    if size(bpop)[0] == 1 and size(bpop)[1] == poped_db["parameters"]["nbpop"]:
        bpop_descr = zeros(poped_db["parameters"]["nbpop"], 3)

        # reorder named values
        bpop = reorder_vec(bpop, fg_names)

        bpop_descr.set_data()[:,1] = bpop
        bpop_descr.set_data()[:,0] = 0  # point values
        bpop_descr.set_data()[:,2] = 0  # variance
        bpop_descr.set_colnam(bpop.get_datanam())
        # bpop_descr.columns.values = bpop.keys()
        bpop = bpop_descr

    # we have just the diagonal parameter values
    if size(sigma)[0] == 1 and (type(sigma) is not np.ndarray or type(sigma) is not matrix):
        sigma_tmp = matrix(np.diag(sigma, size(sigma)[1], size(sigma)[1]))
        sigma_tmp.set_colnam(sigma.keys())
        # sigma_tmp.columns.values = sigma.keys()
        sigma = sigma_tmp

    covd = poped_choose(covd, zeros(
        1, poped_db["parameters"]["NumRanEff"])*(poped_db["parameters"]["NumRanEff"]-1)/2, 0)
    poped_db["parameters"]["covd"] = covd

    tmp = ones(1, len(covd))
    if tmp is not None:
        for i in range(0, len(covd)):
            if covd[i] == 0:
                tmp[i] = 0
    #tmp[covd==0] = 0
    poped_db["parameters"]["notfixed_covd"] = poped_choose(notfixed_covd, tmp, 0)

    # ==================================
    # Sample the individual eta's for FOCE and FOCEI
    # ==================================
    if poped_db["settings"]["iApproximationMethod"] != 0 and poped_db["settings"]["iApproximationMethod"] != 3:

        iMaxCorrIndNeeded = 100
        bzeros = zeros(poped_db["parameters"]["NumRanEff"], 1)
        bones = matrix(np.ones(int(poped_db["parameters"]["NumRanEff"])), (int(poped_db["parameters"]["NumRanEff"]), 1))
        bocczeros = zeros(poped_db["parameters"]["NumDocc"], 1)
        boccones = matrix(np.ones(int(poped_db["parameters"]["NumDocc"])), (int(poped_db["parameters"]["NumDocc"]), 1))

        poped_db["parameters"]["b_global"] = zeros(poped_db["parameters"]["NumRanEff"], max(
            poped_db["settings"]["iFOCENumInd"], iMaxCorrIndNeeded))

        fulld = getfulld(matrix(d.get_data()[:,1]), poped_db["parameters"]["covd"])
        fulldocc = getfulld(matrix(docc.get_data()[:,1]), poped_db["parameters"]["covdocc"])

        poped_db["parameters"]["bocc_global"] = cell(
            poped_db["settings"]["iFOCENumInd"], 1)

        if poped_db["settings"]["d_switch"] is True:
            poped_db["parameters"]["b_global"] = matrix(np.transpose(np.random.multivariate_normal(
                max(poped_db["settings"]["iFOCENumInd"], iMaxCorrIndNeeded), sigma=fulld)))
            for i in range(0, poped_db["settings"]["iFOCENumInd"]):
                poped_db["parameters"]["bocc_global"][i] = zeros(
                    size(docc)[0], poped_db["parameters"]["NumOcc"])
                if poped_db["parameters"]["NumOcc"] != 0:
                    poped_db["parameters"]["bocc_global"][i] = matrix(np.transpose(np.random.multivariate_normal(
                        poped_db["parameters"]["NumOcc"], sigma=fulldocc)))

        else:
            d_dist = pargen(d, poped_db["model"]["user_distribution_pointer"], max(
                poped_db["settings"]["iFOCENumInd"], iMaxCorrIndNeeded), poped_db["settings"]["bLHS"], zeros(1, 0), poped_db)
            docc_dist = pargen(docc, poped_db["model"]["user_distribution_pointer"], poped_db["settings"]
                               ["iFOCENumInd"], poped_db["settings"]["bLHS"], zeros(1, 0), poped_db)

            if len(d_dist) != 0:
                for i in range(0, max(poped_db["settings"]["iFOCENumInd"], iMaxCorrIndNeeded)):
                    poped_db["parameters"]["b_global"][:, i] = matrix(np.transpose(np.random.multivariate_normal(
                        1, sigma=getfulld(d_dist[i, ], poped_db["parameters"]["covd"]))))

            if len(docc_dist) != 0:
                for i in range(0, poped_db["settings"]["iFOCENumInd"]):
                    poped_db["parameters"]["bocc_global"][i] = matrix(np.transpose(np.random.multivariate_normal(
                        poped_db["parameters"]["NumOcc"], sigma=getfulld(docc_dist[i, ], poped_db["parameters"]["covdocc"]))))
    else:
        poped_db["parameters"]["bocc_global"] = zeros(poped_db["parameters"]["NumRanEff"], 1)
        poped_db["parameters"]["bocc_global"] = cell(1, 1)
        poped_db["parameters"]["bocc_global"][1] = zeros(size(docc)[0], poped_db["parameters"]["NumOcc"])
        poped_db["settings"]["iFOCENumInd"] = 1

    poped_db["settings"]["modtit"] = modtit
    poped_db["settings"]["exptit"] = ('%s_exp["mat"]', modtit)  # experiment settings title
    poped_db["settings"]["opttit"] = ('%s_opt["mat"]', modtit)  # optimization settings title
    poped_db["settings"]["bShowGraphs"] = bShowGraphs

    poped_db["settings"]["use_logfile"] = use_logfile
    poped_db["settings"]["output_file"] = output_file
    poped_db["settings"]["output_function_file"] = output_function_file

    poped_db["settings"]["optsw"] = optsw

    line_opta = poped_choose(line_opta, ones(1, size(poped_db["design"]["a"])[1]), 0)
    if test_mat_size(np.array([1, size(poped_db["design"]["a"])[0]]), np.array(line_opta), "line_opta"):
        poped_db["settings"]["line_opta"] = line_opta

    line_optx = poped_choose(line_optx, ones(1, size(poped_db["design"]["x"])[1]), 0)
    if test_mat_size(np.array([1, size(poped_db["design"]["x"])[0]]), np.array(line_opta), "line_optx"):
        poped_db["settings"]["line_optx"] = line_optx


# poped_db$design = popedInput$design
    poped_db["parameters"]["bpop"] = bpop
    poped_db["parameters"]["d"] = d
    poped_db["parameters"]["covd"] = covd
    poped_db["parameters"]["sigma"] = sigma
    poped_db["parameters"]["docc"] = docc
    poped_db["parameters"]["covdocc"] = covdocc

# poped_db["design"]G = design_space$G_xt ## should only be in design_space
# poped_db["design"]Ga = design_space$G_a ## should only be in design_space
# poped_db["design"]Gx = design_space$G_x ## should only be in design_space

# poped_db["design"]maxgroupsize = design_space$maxgroupsize ## should only be in design_space
# poped_db["design"]mingroupsize = design_space$mingroupsize ## should only be in design_space
# poped_db["design"]maxtotgroupsize = design_space$maxtotgroupsize ## should only be in design_space
# poped_db["design"]mintotgroupsize = design_space$mintotgroupsize ## should only be in design_space
# poped_db["design"]maxxt = design_space$maxxt ## should only be in design_space
# poped_db["design"]minxt = design_space$minxt ## should only be in design_space
# poped_db["design"]discrete_x = design_space$discrete_x ## should only be in design_space
# poped_db["design"]maxa = design_space$maxa ## should only be in design_space
# poped_db["design"]mina = design_space$mina ## should only be in design_space

    poped_db["settings"]["m1_switch"] = m1_switch
    poped_db["settings"]["m2_switch"] = m2_switch
    poped_db["settings"]["hle_switch"] = hle_switch
    poped_db["settings"]["gradff_switch"] = gradff_switch
    poped_db["settings"]["gradfg_switch"] = gradfg_switch
    poped_db["settings"]["grad_all_switch"] = grad_all_switch

    poped_db["settings"]["prior_fim"] = prior_fim

    poped_db["parameters"]["notfixed_sigma"] = notfixed_sigma
    # poped_db["parameters"]sigma = sigma
    # poped_db["parameters"]docc = docc

    poped_db["parameters"]["ds_index"] = ds_index

    # create ds_index if not already done
    if poped_db["parameters"]["ds_index"] is not None:
        unfixed_params = get_unfixed_params(poped_db)
        poped_db["parameters"]["ds_index"] = matrix(np.transpose(
            np.repeat(0, unfixed_params["all"])))
        poped_db["parameters"]["ds_index"][(
            unfixed_params["bpop"].get_size() + 1):poped_db["parameters"]["ds_index"].get_size()] = 1
    else:
        if type(poped_db["parameters"]["ds_index"]) is not matrix:
            poped_db["parameters"]["ds_index"] = matrix(np.array(
                [poped_db["parameters"]["ds_index"], 1, len(poped_db["parameters"]["ds_index"])]))

    poped_db["settings"]["strIterationFileName"] = strIterationFileName
    poped_db["settings"]["user_data"] = user_data
    poped_db["settings"]["bUseSecondOrder"] = False
    poped_db["settings"]["bCalculateEBE"] = False
    poped_db["settings"]["bGreedyGroupOpt"] = bGreedyGroupOpt
    poped_db["settings"]["bEANoReplicates"] = bEANoReplicates

    if len(strRunFile) != 0:
        if strRunFile == " ":
            poped_db["settings"]["run_file_pointer"] = zeros(1, 0)
        else:
            if strRunFile is not None:
                poped_db["settings"]["run_file_pointer"] = strRunFile
            else:
                exec(open(popedInput["strRunFile"]).read())
                returnArgs = fileparts(popedInput["strRunFile"])
                strRunFilePath = list(returnArgs.values())[0]
                strRunFilename = list(returnArgs.values())[1]
                # if (~strcmp(strErrorModelFilePath,''))
                # cd(strErrorModelFilePath);
                # end
                poped_db["settings"]["run_file_pointer"] = strRunFilename

#poped_db = convert_popedInput(popedInput,...)

    poped_db["settings"]["Engine"] = {"Type": 1, "Version": poped_version["version"]}

    poped_db = convert_variables(poped_db)  # need to transform here

    param_val = get_all_params(poped_db)
    tmp_names = param_val.keys()
    eval('%s.val = param.val["%s"]' % (tmp_names, tmp_names))
    """ 
	#tools/spped_check.R
	d_val = d_val # for package check
	covd_val = covd_val
	docc_val = docc_val
	covdocc_val = covdocc_val
	bpop_val = bpop_val
	d_full = getfulld(d_val,covd_val)
	docc_full = getfulld(docc_val,covdocc_val)
	sigma_full = poped_db["parameters"]["sigma"]
	
	poped_db["parameters"]["param_pt_val"]["bpop"] = bpop_val
	poped_db["parameters"]["param_pt_val"]["d"] = d_full
	poped_db["parameters"]["param_pt_val"]["docc"] = docc_full
	poped_db["parameters"]["param_pt_val"]["sigma"] = sigma_full
    
	"""
    #     downsize.list = downsizing_general_design(poped_db)
    #     tmp.names = names(downsize.list)
    #     model_switch = c()
    #     ni = c()
    #     xt = c()
    #     x = c()
    #     a = c()
    #     eval(parse(text=paste(tmp.names,"=","downsize.list$",tmp.names,sep="")))
    #     poped_db$downsized.design$model_switch = model_switch
    #     poped_db$downsized.design$ni = ni
    #     poped_db$downsized.design$xt = xt
    #     poped_db$downsized.design$x = x
    #     poped_db$downsized.design$a = a
    #     poped_db$downsized.design$groupsize = poped_db["design"]groupsize

    retargs = fileparts(poped_db["settings"]["output_file"])
    poped_db["settings"]["strOutputFilePath"] = list(retargs.values())[0]
    poped_db["settings"]["strOutputFileName"] = list(retargs.values())[1]
    poped_db["settings"]["strOutputFileExtension"] = list(retargs.values())[2]

    return poped_db

def somestring(**kwargs):
    return ", ".join(f"{key}={value}" for key, value in kwargs.items())

