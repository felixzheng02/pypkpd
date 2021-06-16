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
 Row vector of optimization tasks (1=TRUE,0=FALSE) in the following order: 
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
@param bUseGrouped_xt Use grouped time points (1=TRUE, 0=FALSE).
@param G_xt Matrix defining the grouping of sample points. Matching integers mean that the points are matched.
@param bUseGrouped_a Use grouped covariates (1=TRUE, 0=FALSE)
@param G_a Matrix defining the grouping of covariates. Matching integers mean that the points are matched.
@param bUseGrouped_x Use grouped discrete design variables (1=TRUE, 0=FALSE).
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
\item 7 = Inverse of the sum of the expected parameter RSE: 1/sum(get_rse(FIM,poped.db,use_percent=FALSE))
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
That is, from your full IIV matrix  \code{covd <-  IIV[lower.tri(IIV)]}. 
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
Use random search (1=TRUE, 0=FALSE)
@param bUseStochasticGradient Use Stochastic Gradient search (1=TRUE, 0=FALSE) 
@param bUseLineSearch Use Line search (1=TRUE, 0=FALSE) 
@param bUseExchangeAlgorithm Use Exchange algorithm (1=TRUE, 0=FALSE)        
@param bUseBFGSMinimizer Use BFGS Minimizer (1=TRUE, 0=FALSE) 
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
@param dSeed The seed number used for optimization and sampling -- integer or -1 which creates a random seed \code{as.integer(Sys.time())} or NULL.
@param line_opta Vector for line search on continuous design variables (1=TRUE,0=FALSE)
@param line_optx Vector for line search on discrete design variables (1=TRUE,0=FALSE) 
@param bShowGraphs Use graph output during search
@param use_logfile If a log file should be used (0=FALSE, 1=TRUE)
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
@param iDiffSolverMethod The diff equation solver method, NULL as default.
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


#from am.poped_choose import am.poped_choose
import project.all_modules as am



def create_poped_database(
						popedInput = {},
		):

	# parameters

	## --------------------------
	## ---- Model definition
	## --------------------------

	# -- Filname and path of the model file --
	ff_file = None
	ff_fun = am.poped_choose(popedInput['model']['ff_pointer'], None, 0)
	# -- Filname and path of the g parameter file --
	fg_file = None
	fg_fun = am.poped_choose(popedInput['model']['fg_pointer'], None, 0)
	# -- Filname and path of the error model file --
	fError_file = None
	fError_fun = am.poped_choose(popedInput['model']['ferror_pointer'], None, 0)

	## --------------------------
	## ---- What to optimize
	## --------------------------

	## -- Vector of optimization tasks (1=TRUE,0=FALSE)
	## (Samples per subject, Sampling schedule, Discrete design variable, Continuous design variable, Number of id per group)
	## -- All elements set to zero => only calculate the FIM with current design --
	optsw = am.poped_choose(popedInput['settings']['optsw'], am.np.array([[0,0,0,0,0]]), 0)

	## --------------------------
	## ---- Initial Design
	## --------------------------

	## -- Matrix defining the initial sampling schedule --
	xt = am.poped_choose(popedInput['design'][['xt']],"'xt' needs to be defined", 1)
	## -- Number of groups/individuals --
	#size(,1) find number of rows
	#m=am.poped_choose(popedInput[['m']], size(xt,1)),
	m = am.poped_choose(popedInput['design'][['m']], None, 0)
	## -- Matrix defining the initial discrete values --
	#x=am.poped_choose(popedInput['design'][['x']], zeros(m,0)),
	x = am.poped_choose(popedInput['design'][['x']], None, 0)
	## -- Number of discrete design variables --
	#size(,1) find number of cols
	#nx=am.poped_choose(popedInput['nx'], size(x,2)),
	nx = am.poped_choose(popedInput['design']['nx'], None, 0)
	## -- Vector defining the initial covariate values --
	#a=poped.choose(popedInput$design[["a"]],zeros(m,0)),
	a = am.poped_choose(popedInput['design'][['a']], None, 0)
	## number of continuous design variables that are not time (e.g. continuous covariates)
	#na=am.poped_choose(popedInput['na'],size(a,2)),
	#na=am.poped_choose(popedInput['na'],None),
	## -- Vector defining the size of the different groups (num individuals in each group) --
	groupsize = am.poped_choose(popedInput['design']['groupsize'], "'groupsize' needs to be defined", 1)
	## -- Vector defining the number of samples for each group --
	#ni=poped.choose(popedInput$design$ni,matrix(size(xt,2),m,1)),
	ni = am.poped_choose(popedInput['design']['ni'], None, 0)
	## -- Vector defining which response a certain sampling time belongs to --
	#model_switch=poped.choose(popedInput$design$model_switch,ones(size(xt,1),size(xt,2))),
	model_switch = am.poped_choose(popedInput['design']['model_switch'], None, 0)

	## --------------------------
	## ---- Design space
	## --------------------------

	## -- Max number of samples per group/individual --
	maxni = am.poped_choose(popedInput['design_space']['maxni'], None, 0)
	## -- Min number of samples per group/individual --
	minni = am.poped_choose(popedInput['design_space']['minni'], None, 0)
	maxtotni = am.poped_choose(popedInput['design_space']['maxtotni'], None, 0)
	mintotni = am.poped_choose(popedInput['design_space']['mintotni'], None, 0)
	## -- Vector defining the max size of the different groups (max num individuals in each group) --
	maxgroupsize = am.poped_choose(popedInput['design_space']['maxgroupsize'], None, 0)
	## -- Vector defining the min size of the different groups (min num individuals in each group) --
	#mingroupsize=poped.choose(popedInput$design$mingroupsize,ones(m,1)),
	mingroupsize = am.poped_choose(popedInput['design_space']['mingroupsize'], None, 0)
	## -- The total maximal groupsize over all groups--
	maxtotgroupsize = am.poped_choose(popedInput['design_space']['maxtotgroupsize'], None, 0)
	## -- The total minimal groupsize over all groups--
	mintotgroupsize = am.poped_choose(popedInput['design_space']['mintotgroupsize'], None, 0)
	## -- Matrix defining the max value for each sample --
	maxxt = am.poped_choose(popedInput['design_space']['maxxt'], None, 0)
	## -- Matrix defining the min value for each sample --
	minxt = am.poped_choose(popedInput['design_space']['minxt'], None, 0)
	discrete_xt = am.poped_choose(popedInput['design_space']['xt_space'], None, 0)
	## -- Cell defining the discrete variables --
	#discrete_x=poped.choose(popedInput$design$discrete_x,cell(m,nx)),
	discrete_x = am.poped_choose(popedInput['design_space']['discrete_x'], None, 0)
	## -- Vector defining the max value for each covariate --
	maxa = am.poped_choose(popedInput['design_space']['maxa'], None, 0)
	## -- Vector defining the min value for each covariate --
	mina = am.poped_choose(popedInput['design_space']['mina'], None, 0)
	discrete_a = am.poped_choose(popedInput['design_space']['a_space'], None, 0)
	## -- Use grouped time points (1=TRUE, 0=FALSE) --
	bUseGrouped_xt = am.poped_choose(popedInput['design_space']['bUseGrouped_xt'], False, 0)
	## -- Matrix defining the grouping of sample points --
	G_xt = am.poped_choose(popedInput['design_space']['G_xt'], None, 0)
	## -- Use grouped covariates (1=TRUE, 0=FALSE) --
	bUseGrouped_a = am.poped_choose(popedInput['design_space']['bUseGrouped_a'], False, 0)
	## -- Matrix defining the grouping of covariates --
	G_a = am.poped_choose(popedInput['design_space']['G_a'], None, 0)
	## -- Use grouped discrete design variables (1=TRUE, 0=FALSE) --
	bUseGrouped_x = am.poped_choose(popedInput['design_space']['bUseGrouped_x'], False, 0)
	## -- Matrix defining the grouping of discrete design variables --
	G_x = am.poped_choose(popedInput['design_space'][["G_x"]], None, 0)

	## --------------------------
	## ---- FIM calculation
	## --------------------------

	## -- Fisher Information Matrix type
	## (0=Full FIM,
	##  1=Reduced FIM,
	##  2=weighted models,
	##  3=Loc models,
	##  4=reduced FIM with derivative of SD of sigma as pfim,
	##  5=FULL FIM parameterized with A,B,C matrices & derivative of variance,
	##  6=Calculate one model switch at a time, good for large matrices,
	##  7=Reduced FIM parameterized with A,B,C matrices & derivative of variance) --
	iFIMCalculationType = am.poped_choose(popedInput['settings']['iFIMCalculationType'], 1, 0)
	## -- Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI --
	iApproximationMethod = am.poped_choose(popedInput['settings']['iApproximationMethod'], 0, 0)
	## -- Num individuals in each step of FOCE --
	iFOCENumInd = am.poped_choose(popedInput['settings']['iFOCENumInd'], 1000, 0)
	## -- The prior FIM (added to calculated FIM) --
	prior_fim = am.poped_choose(popedInput['settings']['prior_fim'], am.np.array([0]), 0)
	## -- Filname and path for the Autocorrelation function, empty string means no autocorrelation --
	strAutoCorrelationFile = am.poped_choose(popedInput['model']['auto_pointer'], "", 0)

	## --------------------------
	## ---- Criterion specification
	## --------------------------

	## -- D-family design (1) or ED-family design (0) (with or without parameter uncertainty) --
	d_switch = am.poped_choose(popedInput['settings']['d_switch'], 1, 0)
	## -- OFV calculation type for FIM (1=Determinant of FIM,4=log determinant of FIM,6=determinant of interesting part of FIM (Ds)) --
	ofv_calc_type = am.poped_choose(popedInput['settings']['ofv_calc_type'], 4, 0)
	## -- Ds_index, set index to 1 if a parameter is uninteresting, otherwise 0.
	## size=(1,num unfixed parameters). First unfixed bpop, then unfixed d, then unfixed docc and last unfixed sigma --
	## default is the fixed effects being important
	ds_index = popedInput['parameters']['ds_index']
	## -- Penalty function, empty string means no penalty.  User defined criterion --
	strEDPenaltyFile = am.poped_choose(popedInput['settings']['strEDPenaltyFile'], "", 0)
	ofv_fun = am.poped_choose(popedInput['settings']['ofv_fun'], None, 0)



	## --------------------------
	## ---- E-family Criterion options
	## --------------------------
	## -- ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
	iEDCalculationType = am.poped_choose(popedInput['settings']['iEDCalculationType'], 0, 0)
	## -- Sample size for E-family sampling --
	ED_samp_size = am.poped_choose(popedInput['settings']['ED_samp_size'], 45, 0)
	## -- How to sample from distributions in E-family calculations. 0=Random Sampling, 1=LatinHyperCube --
	bLHS = am.poped_choose(popedInput['settings']['bLHS'], 1, 0)
	## -- Filname and path for user defined distributions for E-family designs --
	strUserDistributionFile = am.poped_choose(popedInput['model']['user_distribution_pointer'], "", 0)

	## --------------------------
	## ---- Model parameters
	## --------------------------

	## -- Number of typical values --
	nbpop = popedInput['parameters']['nbpop']
	## -- Number of IIV parameters --
	NumRanEff = popedInput['parameters']['NumRanEff']
	## -- Number of IOV variance parameters --
	NumDocc = popedInput['parameters']['NumDocc']
	## -- Number of occassions --
	NumOcc = popedInput['parameters']['NumOcc']
	## -- The length of the g parameter vector --
	#ng=popedInput["parameters"]ng,

	## -- Matrix defining the fixed effects, per row (row number = parameter_number),
	## the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
	## 3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal).
	## The second column defines the mean.
	## The third column defines the variance of the distribution.
	# can also just supply the parameter values as a c()
	bpop = am.poped_choose(popedInput['parameters']['bpop'], 'bpop must be defined', 1)
	## -- Matrix defining the diagonals of the IIV (same logic as for the fixed effects) --
	# can also just supply the parameter values as a c()
	d = am.poped_choose(popedInput['parameters']['d'], None, 0)	## -- vector defining the row major lower triangle of the covariances of the IIV variances --
	# set to zero if not defined
	covd = popedInput['parameters']['covd']
	## -- Matrix defining the variances of the residual variability terms --
	## REQUIRED! No defaults given.
	# can also just supply the diagonal values as a c()
	sigma: am.np.ndarray = popedInput['parameters']['sigma']
	## -- Matrix defining the IOV, the IOV variances and the IOV distribution --
	docc = am.poped_choose(popedInput['parameters']['docc'], am.np.array([0, 0, 0]), 0)
	## -- Matrix defining the covariance of the IOV --
	covdocc = am.poped_choose(popedInput['parameters']['covdocc'], am.zeros(1, len(docc[:, 1])*(len(docc[:, 2])-1)/2), 0)

	## --------------------------
	## ---- Model parameters fixed or not
	## --------------------------
	## -- Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed) --
	notfixed_bpop = popedInput['parameters']['notfixed_bpop']
	## -- Vector defining if a IIV is fixed or not (1=not fixed, 0=fixed) --
	notfixed_d = popedInput['parameters']['notfixed_d']
	## -- Vector defining if a covariance IIV is fixed or not (1=not fixed, 0=fixed) --
	notfixed_covd = popedInput['parameters']['notfixed_covd']
	## -- Vector defining if an IOV variance is fixed or not (1=not fixed, 0=fixed) --
	notfixed_docc = popedInput['["parameters"]']['notfixed_docc']
	## -- Vector row major order for lower triangular matrix defining if a covariance IOV is fixed or not (1=not fixed, 0=fixed) --
	notfixed_covdocc = am.poped_choose(popedInput['parameters']['notfixed_covdocc'], am.zeros(1, len(covdocc)), 0)
	## -- Vector defining if a residual error parameter is fixed or not (1=not fixed, 0=fixed) --
	notfixed_sigma = am.poped_choose(popedInput['parameters']['notfixed_sigma'], am.np.ones(sigma.shape[1]), 0)
	## -- Vector defining if a covariance residual error parameter is fixed or not (1=not fixed, 0=fixed) --
	## default is fixed
	notfixed_covsigma = am.poped_choose(popedInput['parameters']['notfixed_covsigma'], am.zeros(1,len(notfixed_sigma)*(len(notfixed_sigma)-1)/2), 0)


	## --------------------------
	## ---- Optimization algorithm choices
	## --------------------------

	## -- Use random search (1=TRUE, 0=FALSE) --
	bUseRandomSearch = am.poped_choose(popedInput['settings']['bUseRandomSearch'], True, 0)
	## -- Use Stochastic Gradient search (1=TRUE, 0=FALSE) --
	bUseStochasticGradient = am.poped_choose(popedInput['settings']['bUseStochasticGradient'], True, 0)
	## -- Use Line search (1=TRUE, 0=FALSE) --
	bUseLineSearch = am.poped_choose(popedInput['settings']['bUseLineSearch'], True, 0)
	## -- Use Exchange algorithm (1=TRUE, 0=FALSE) --
	bUseExchangeAlgorithm = am.poped_choose(popedInput['settings']['bUseExchangeAlgorithm'], False, 0)
	## -- Use BFGS Minimizer (1=TRUE, 0=FALSE) --
	bUseBFGSMinimizer = am.poped_choose(popedInput['settings']['bUseBFGSMinimizer'], False, 0)
	## -- Exchange Algorithm Criteria, 1 = Modified, 2 = Fedorov --
	EACriteria = am.poped_choose(popedInput['settings']['EACriteria'], 1, 0)
	## -- Filename and path for a run file that is used instead of the regular PopED call --
	strRunFile = am.poped_choose(popedInput['settings']['run_file_pointer'], "", 0)

	## --------------------------
	## ---- Labeling and file names
	## --------------------------

	## -- The current PopED version --
# ！！！没写 packageVersion("PopED")
	poped_version = am.poped_choose(popedInput['settings']['poped_version'], "0.0.2", 0)
	## -- The model title --
	modtit = am.poped_choose(popedInput['settings']['modtit'],'PopED model', 0)
	## -- Filname and path of the output file during search --
	output_file = am.poped_choose(popedInput['settings']['output_file'], "PopED_output_summary", 0)
	## -- Filname suffix of the result function file --
	output_function_file = am.poped_choose(popedInput['settings']['output_function_file'], "PopED_output_", 0)
	## -- Filename and path for storage of current optimal design --
	strIterationFileName = am.poped_choose(popedInput['settings']['strIterationFileName'], "PopED_current.R", 0)


	## --------------------------
	## ---- Misc options
	## --------------------------
	## -- User defined data structure that, for example could be used to send in data to the model --
	user_data = am.poped_choose(popedInput['settings']['user_data'], am.cell(0,0), 0)
	## -- Value to interpret as zero in design --
	ourzero = am.poped_choose(popedInput['settings']['ourzero'], 1e-5, 0)
	#ourzero=poped.choose(popedInput$ourzero,0),
	## -- The seed number used for optimization and sampling -- integer or -1 which creates a random seed
	dSeed = am.poped_choose(popedInput['settings']['dSeed'], None, 0)
	## -- Vector for line search on continuous design variables (1=TRUE,0=FALSE) --
	line_opta = am.poped_choose(popedInput['settings']['line_opta'], None, 0)
	## -- Vector for line search on discrete design variables (1=TRUE,0=FALSE) --
	line_optx = am.poped_choose(popedInput['settings']['line_optx'], None, 0) #matrix(0,0,1)
	## -- Use graph output during search --
	bShowGraphs = am.poped_choose(popedInput['settings']['bShowGraphs'], False, 0)
	## -- If a log file should be used (0=FALSE, 1=TRUE) --
	use_logfile = am.poped_choose(popedInput['settings']['use_logfile'], False, 0)
	## -- Method used to calculate M1
	## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
	m1_switch = am.poped_choose(popedInput['settings']['m1_switch'], 1, 0)
	## -- Method used to calculate M2
	## (0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
	m2_switch = am.poped_choose(popedInput['settings']['m2_switch'], 1, 0)
	## -- Method used to calculate linearization of residual error
	## (0=Complex difference, 1=Central difference, 30=Automatic differentiation) --
	hle_switch = am.poped_choose(popedInput['settings']['hle_switch'], 1, 0)
	## -- Method used to calculate the gradient of the model
	## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
	gradff_switch = am.poped_choose(popedInput['settings']['gradff_switch'], 1, 0)
	## -- Method used to calculate the gradient of the parameter vector g
	## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
	gradfg_switch = am.poped_choose(popedInput['settings']['gradfg_switch'], 1, 0)
	## -- Method used to calculate all the gradients
	## (0=Complex difference, 1=Central difference) --
	grad_all_switch = am.poped_choose(popedInput['settings']['grad_all_switch'], 1, 0)
	## -- Number of iterations in random search between screen output --
	rsit_output = am.poped_choose(popedInput['settings']['rsit_output'], 5, 0)
	## -- Number of iterations in stochastic gradient search between screen output --
	sgit_output = am.poped_choose(popedInput['settings']['sgit_output'], 1, 0)
	## -- Step length of derivative of linearized model w.r.t. typical values --
	hm1 = am.poped_choose(popedInput['settings'][["hm1"]], 0.00001, 0)
	## -- Step length of derivative of model w.r.t. g --
	hlf = am.poped_choose(popedInput['settings'][["hlf"]], 0.00001, 0)
	## -- Step length of derivative of g w.r.t. b --
	hlg = am.poped_choose(popedInput['settings'][["hlg"]], 0.00001, 0)
	## -- Step length of derivative of variance w.r.t. typical values --
	hm2 = am.poped_choose(popedInput['settings'][["hm2"]], 0.00001, 0)
	## -- Step length of derivative of OFV w.r.t. time --
	hgd = am.poped_choose(popedInput['settings'][["hgd"]], 0.00001, 0)
	## -- Step length of derivative of model w.r.t. sigma --
	hle = am.poped_choose(popedInput['settings'][["hle"]], 0.00001, 0)
	## -- The absolute tolerance for the diff equation solver --
	AbsTol = am.poped_choose(popedInput['settings']['AbsTol'], 0.000001, 0)
	## -- The relative tolerance for the diff equation solver --
	RelTol = am.poped_choose(popedInput['settings']['RelTol'], 0.000001, 0)
	## -- The diff equation solver method, 0, no other option --
	iDiffSolverMethod = am.poped_choose(popedInput['settings']['iDiffSolverMethod'], None, 0)
	## -- If the differential equation results should be stored in memory (1) or not (0) --
	bUseMemorySolver = am.poped_choose(popedInput['settings']['bUseMemorySolver'], False, 0)
	## -- Number of Random search iterations --
	rsit = am.poped_choose(popedInput['settings'][["rsit"]], 300, 0)
	## -- Number of Stochastic gradient search iterations --
	sgit = am.poped_choose(popedInput['settings'][["sgit"]], 150, 0)
	## -- Number of Random search iterations with discrete optimization --
	intrsit = am.poped_choose(popedInput['settings']['intrsit'], 250, 0)
	## -- Number of Stochastic Gradient search iterations with discrete optimization --
	intsgit = am.poped_choose(popedInput['settings']['intsgit'], 50, 0)
	## -- Iterations until adaptive narrowing in random search --
	maxrsnullit = am.poped_choose(popedInput['settings']['maxrsnullit'], 50, 0)
	## -- Stoachstic Gradient convergence value,
	## (difference in OFV for D-optimal, difference in gradient for ED-optimal) --
	convergence_eps = am.poped_choose(popedInput['settings']['convergence_eps'], 1e-08, 0)
	## -- Random search locality factor for sample times --
	rslxt = am.poped_choose(popedInput['settings']['rslxt'], 10, 0)
	## -- Random search locality factor for covariates --
	rsla = am.poped_choose(popedInput['settings']['rsla'], 10, 0)
	## -- Stochastic Gradient search first step factor for sample times --
	cfaxt = am.poped_choose(popedInput['settings']['cfaxt'], 0.001, 0)
	## -- Stochastic Gradient search first step factor for covariates --
	cfaa = am.poped_choose(popedInput['settings']['cfaa'], 0.001, 0)
	## -- Use greedy algorithm for group assignment optimization --
	bGreedyGroupOpt = am.poped_choose(popedInput['settings']['bGreedyGroupOpt'], False, 0)
	## -- Exchange Algorithm StepSize --
	EAStepSize = am.poped_choose(popedInput['settings']['EAStepSize'], 0.01, 0)
	## -- Exchange Algorithm NumPoints --
	EANumPoints = am.poped_choose(popedInput['settings']['EANumPoints'], False, 0)
	## -- Exchange Algorithm Convergence Limit/Criteria --
	EAConvergenceCriteria = am.poped_choose(popedInput['settings']['EAConvergenceCriteria'], 1e-20, 0)
	## -- Avoid replicate samples when using Exchange Algorithm --
	bEANoReplicates = am.poped_choose(popedInput['settings']['bEANoReplicates'], False, 0)
	## -- BFGS Minimizer Convergence Criteria Minimum Step --
	BFGSConvergenceCriteriaMinStep = None,
	#poped.choose(popedInput["settings"]BFGSConvergenceCriteriaMinStep,1e-08),
	## -- BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance --
	BFGSProjectedGradientTol = am.poped_choose(popedInput['settings']['BFGSProjectedGradientTol'], 0.0001, 0)
	## -- BFGS Minimizer Line Search Tolerance f --
	BFGSTolerancef = am.poped_choose(popedInput['settings']['BFGSTolerancef'], 0.001, 0)
	## -- BFGS Minimizer Line Search Tolerance g --
	BFGSToleranceg = am.poped_choose(popedInput['settings']['BFGSToleranceg'], 0.9, 0)
	## -- BFGS Minimizer Line Search Tolerance x --
	BFGSTolerancex = am.poped_choose(popedInput['settings']['BFGSTolerancex'], 0.1, 0)
	## -- Number of iterations in ED-optimal design to calculate convergence criteria --
	ED_diff_it = am.poped_choose(popedInput['settings']['ED_diff_it'], 30, 0)
	## -- ED-optimal design convergence criteria in percent --
	ED_diff_percent = am.poped_choose(popedInput['settings']['ED_diff_percent'], 10, 0)
	## -- Number of grid points in the line search --
	line_search_it = am.poped_choose(popedInput['settings']['ls_step_size'], 50, 0)
	## -- Number of iterations of full Random search and full Stochastic Gradient if line search is not used --
	Doptim_iter = am.poped_choose(popedInput['settings']['iNumSearchIterationsIfNotLineSearch'], 1, 0)

	## --------------------------
	## -- Parallel options for PopED -- --
	## --------------------------
	#     ## -- Compile option for PopED
	#     ## -1 = No compilation,
	#     ## 0 or 3 = Full compilation,
	#     ## 1 or 4 = Only using MCC (shared lib),
	#     ## 2 or 5 = Only MPI,
	#     ## Option 0,1,2 runs PopED and option 3,4,5 stops after compilation --
	iCompileOption = am.poped_choose(popedInput['settings']['parallel']['iCompileOption'], -1, 0)
	## -- Parallel method to use (0 = Matlab PCT, 1 = MPI) --
	iUseParallelMethod = am.poped_choose(popedInput['settings']['parallel']['iUseParallelMethod'], 1, 0)
	## -- Additional dependencies used in MCC compilation (mat-files), if several space separated --
	MCC_Dep = None,
	#poped.choose(popedInput["settings"]parallel$strAdditionalMCCCompilerDependencies, ''),
	## -- Compilation output executable name --
	strExecuteName = am.poped_choose(popedInput['settings']['parallel']['strExecuteName'], 'calc_fim.exe', 0)
	## -- Number of processes to use when running in parallel (e.g. 3 = 2 workers, 1 job manager) --
	iNumProcesses = am.poped_choose(popedInput['settings']['parallel']['iNumProcesses'], 2, 0)
	## -- Number of design evaluations that should be evaluated in each process before getting new work from job manager --
	iNumChunkDesignEvals = am.poped_choose(popedInput['settings']['parallel']['iNumChunkDesignEvals'], -2, 0)
	## -- The prefix of the input mat file to communicate with the executable --
	#strMatFileInputPrefix = poped.choose(
	#   popedInput["settings"]parallel$strMatFileInputPrefix,
	#   'parallel_input'),
	## -- The prefix of the output mat file to communicate with the executable --
	Mat_Out_Pre = am.poped_choose(popedInput['settings']['parallel']['strMatFileOutputPrefix'], 'parallel_output', 0)
	## -- Extra options send to e$g. the MPI executable or a batch script, see execute_parallel$m for more information and options --
	strExtraRunOptions = am.poped_choose(popedInput['settings']['parallel']['strExtraRunOptions'], '', 0)
	## -- Polling time to check if the parallel execution is finished --
	dPollResultTime = am.poped_choose(popedInput['settings']['parallel']['dPollResultTime'], 0.1, 0)
	## -- The file containing the popedInput structure that should be used to evaluate the designs --
	strFunctionInputName = am.poped_choose(popedInput['settings']['parallel']['strFunctionInputName'],'function_input', 0)
	## -- If the random search is going to be executed in parallel --
	bParallelRS = am.poped_choose(popedInput['settings']['parallel']['bParallelRS'], False, 0)
	## -- If the stochastic gradient search is going to be executed in parallel --
	bParallelSG = am.poped_choose(popedInput['settings']['parallel']['bParallelSG'], False, 0)
	## -- If the modified exchange algorithm is going to be executed in parallel --
	bParallelMFEA = am.poped_choose(popedInput['settings']['parallel']['bParallelMFEA'], False, 0)
	## -- If the line search is going to be executed in parallel --
	bParallelLS = am.poped_choose(popedInput['settings']['parallel']['bParallelLS'], False, 0)

#####-------------- main part --------------#####
	poped_db = {}
    # five main headings for database
    #     poped.db <- list(design=NULL,
    #                      design_space=NULL,
    #                      models=NULL,
    #                      parameters=NULL,
    #                      settings=NULL)
    
    #     # update popedInput with options supplied in function
    #     called_args <- match.call()
    #     default_args <- formals()
    #     for(i in names(called_args)[-1]){
    #       if(length(grep("^popedInput$",capture.output(default_args[[i]])))==1) {
    #         eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
    #       }
    #     }
    
    #modifyList(settings, list()$settings)
    
    ## compare to a default input function.
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
		BFGSConvergenceCriteriaMinStep = am.poped_choose(popedInput["settings"]["BFGSConvergenceCriteriaMinStep"], 1e-08, 0)
	
	if MCC_Dep is None:
		MCC_Dep = am.poped_choose(popedInput["settings"]["parallel"]["strAdditionalMCCCompilerDependencies"], "", 0)
	
	poped_db["model"] = {}
	poped_db["model"]["user_distribution_pointer"] = ""
	
	if str(strUserDistributionFile) != "":
		if strUserDistributionFile != None:
			poped_db["model"]["user_distribution_pointer"] = strUserDistributionFile
		else:
			exec(open(strUserDistributionFile).read)
			returnArgs = am.fileparts(strUserDistributionFile)
			strUserDistFilePath = returnArgs[[0]]
			strUserDistFilename = returnArgs[[1]]
			poped_db["model"]["user_distribution_pointer"] = strUserDistFilename
	
	#should be removed
	if any(x.shape[0]==0 & x.shape[1]==0):
		x = None
		G_x = None
		discrete_x = None
	
	#should be removed
	if any(a.shape[0]==0 & a.shape[1]==0):
		a = None
		G_a = None
		mina = None
		maxa = None

	design = am.create_design(xt=xt,
                        groupsize=groupsize,
                        m=m,
                        x=x,
                        a=a,
                        ni=ni,
                        model_switch=model_switch)
		
	design_space = am.create_design_space(design,
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
                                    xt_space = discrete_xt,
									maxa=maxa,   
                                    mina=mina, 
                                    a_space = discrete_a, 
                                    x_space = discrete_x,    
                                    use_grouped_xt=bUseGrouped_xt, 
                                    grouped_xt=G_xt, 
                                    use_grouped_a=bUseGrouped_a,               
                                    grouped_a=G_a,               
                                    use_grouped_x=bUseGrouped_x,               
                                    grouped_x=G_x,                                        
									our_zero=ourzero)

	design = design_space["design"]
	design_space = design_space["design_space"]

	## all of this should be replaced with using the names used in create_design_space as function arguments
	if design_space[["use_grouped_a"]] != None:
		design_space["bUseGrouped_a"] = design_space[["use_grouped_a"]]
    	design_space[["use_grouped_a"]] = None
    
	if design_space[["use_grouped_x"]] != None:
		design_space["bUseGrouped_x"] = design_space[["use_grouped_x"]]
    	design_space[["use_grouped_x"]] = None
    
	if design_space[["use_grouped_xt"]] != None:
		design_space["bUseGrouped_xt"] = design_space[["use_grouped_xt"]]
    	design_space[["use_grouped_xt"]] = None
    
	if design_space[["grouped_a"]] != None:
		design_space["G_a"] = design_space[["grouped_a"]]
    	design_space[["grouped_a"]] = None
    
	if design_space[["grouped_x"]] != None:
		design_space["G_x"] = design_space[["grouped_x"]]
    	design_space[["grouped_x"]] = None
    
	if design_space[["grouped_xt"]] != None:
		design_space["G_xt"] = design_space[["grouped_xt"]]
    	design_space[["grouped_xt"]] = None
    
	if design_space[["x_space"]] != None:
		design_space["discrete_x"] = design_space[["x_space"]]
    	design_space[["x_space"]] = None
    
    #design_space$maxni <- max(design_space$maxni)
    #design_space$minni <- min(design_space$minni)	
	
	#should be removed
	if design[["x"]] is None:
		design["x"] = am.zeros(design["m"], 0)
		design_space["G_x"] = design["x"]
		design_space["bUseGrouped_x"] = False
		design_space["discrete_x"] = am.cell(design["m"],0)
	
	#should be removed
	if design[["a"]] is None:
		design["a"] = am.zeros(design["m"], 0)
		design_space["G_a"] = design["a"]
		design_space["bUseGrouped_a"] = False
		design_space["mina"] = design["a"]
		design_space["maxa"] = design["a"]

	poped_db["design"] = design
	poped_db["design_space"] = design_space

	#poped.db$m = poped.db$design$m  # should be removed only in design
    #poped.db$nx = poped.choose(nx,size(design$x,2)) # should be removed, not needed or in design
    #poped.db$na = poped.choose(na,size(design$a,2)) # should be removed, not needed or in design
                
	poped_db["settings"]["bLHS"] = bLHS
    
    #poped.db$discrete_x = design_space$discrete_x # should be removed only in design_space
    
    #poped.db$maxni=max(design_space$maxni) # should be only in design_space and called maxmaxni if needed
    #poped.db$minni=min(design_space$minni) # should be only in design_space and called minminni if needed
    
    #poped.db$bUseGrouped_xt = design_space$bUseGrouped_xt # should be only in design_space
    #poped.db$bUseGrouped_a  = design_space$bUseGrouped_a # should be only in design_space
    #poped.db$bUseGrouped_x  = design_space$bUseGrouped_x # should be only in design_space
    
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

	poped_db["parameters"]["covdocc"] = covdocc
	poped_db["parameters"]["notfixed_covdocc"] = notfixed_covdocc
	poped_db["parameters"]["notfixed_covsigma"] = notfixed_covsigma

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
	#Temp thing for memory solvers
	poped_db["settings"]["bUseMemorySolver"] = bUseMemorySolver
	poped_db["settings"]["solved_solutions"] = am.cell(0,0)
	poped_db["settings"]["maxtime"] = max(max(poped_db["design_space"]["maxxt"])) + poped_db["settings"]["hgd"]

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
		if fg_fun != None:
			poped_db["model"]["fg_pointer"] = fg_fun
	elif fg_file != None:
		poped_db["model"]["fg_pointer"] = fg_file
	else:
		exec(open(fg_file).read)
		returnArgs =  am.fileparts(fg_file) 
      	strfgModelFilePath = returnArgs[[0]]
      	strfgModelFilename  = returnArgs[[1]]
      	## if (~strcmp(strfgModelFilePath,''))
      	##    cd(strfgModelFilePath);
      	## end
      	poped_db["model"]["fg_pointer"] = strfgModelFilename
	
	poped_db["settings"]["ed_penalty_pointer"] = am.zeros(1,0)
	if str(strEDPenaltyFile != ""):
		if strEDPenaltyFile != None:
			poped_db["settings"]["ed_penalty_pointer"] = strEDPenaltyFile
		else:
			exec(open(popedInput["strEDPenaltyFile"]).read)
			returnArgs =  am.fileparts(popedInput["strEDPenaltyFile"]) 
        strEDPenaltyFilePath = returnArgs[[0]]
        strEDPenaltyFilename = returnArgs[[1]]
        ##     if (~strcmp(strEDPenaltyFilePath,''))
        ##        cd(strEDPenaltyFilePath);
        ##     end
        poped_db["settings"]["ed_penalty_pointer"] = strEDPenaltyFilename


	# if(is.null(ofv_fun) || is.function(ofv_fun)){
    #   poped.db["settings"]ofv_fun = ofv_fun
    # } else {
    #   stop("ofv_fun must be a function or NULL")
    # }

	if ofv_fun is None or callable(ofv_fun):
		poped_db["settings"]["ofv_fun"] = ofv_fun
	else:
		# source explicit file
        # here I assume that function in file has same name as filename minus .txt and pathnames
		if ofv_fun != None:
			exec(open(str(ofv_fun).read))
			poped_db["settings"]["ofv_fun"] = ????eval(parse(text=fileparts(ofv_fun)[["filename"]]))
		else:
			raise Exception("ofv_fun is not a function or NULL, and no file with that name was found")

	
	# if(is.function(ofv_fun)){
    #   poped.db["settings"]ofv_fun = ofv_fun 
    # } else if(exists(ofv_fun)){
    #   poped.db["settings"]ofv_fun = ofv_fun 
    # } else {
    #   source(ofv_fun)
    #   returnArgs <-  fileparts(ofv_fun) 
    #   strffModelFilePath <- returnArgs[[1]]
    #   strffModelFilename  <- returnArgs[[2]]
    #   ## if (~strcmp(strffModelFilePath,''))
    #   ##    cd(strffModelFilePath);
    #   ## end
    #   poped.db["settings"]ofv_fun = strffModelFilename
    # }
 
	poped_db["model"]["auto_pointer"] = am.zeros(1,0)
	#poped_db["model"]["auto_pointer"] = ''
	if strAutoCorrelationFile != "":
		if str(strAutoCorrelationFile) != "":
			if strAutoCorrelationFile != None:
				poped_db["model"]["auto_pointer"] = strAutoCorrelationFile
			else:
				exec(open(popedInput["strAutoCorrelationFile"]).read) 
				returnArgs =  am.fileparts(popedInput["strAutoCorrelationFile"]) 
				strAutoCorrelationFilePath = returnArgs[[0]]
				strAutoCorrelationFilename  = returnArgs[[1]]
          		##     if (~strcmp(strAutoCorrelationFilePath,''))
         		##        cd(strAutoCorrelationFilePath);
          		##     end
				poped_db["model"]["auto_pointer"] = strAutoCorrelationFilename

	if callable(ff_fun):
		poped_db["model"]["ff_pointer"] = ff_fun
	elif type(ff_fun) is str:
		if ff_fun != None:
			poped_db["model"]["ff_pointer"] = ff_fun
	elif ff_file != None:
		poped_db["model"]["ff_pointer"] = ff_file
	else:
		exec(open(popedInput["strAutoCorrelationFile"]).read)
		returnArgs =  am.fileparts(ff_file) 
		strffModelFilePath = returnArgs[[0]]
		strffModelFilename  = returnArgs[[1]]
		## if (~strcmp(strffModelFilePath,''))
        ##    cd(strffModelFilePath);
        ## end
		poped_db["model"]["ff_pointer"] = strffModelFilename

	#Check if there is any sub models defined
	if any(names(popedInput) == "SubModels"))
      i = 1
      while any(names(popedInput$SubModels)==sprintf('ff_file%d',i)):
        source(eval(sprintf('popedInput$SubModels$ff_file%d',i))) ##ok<NASGU> 
        returnArgs <-  fileparts(eval(sprintf('popedInput$SubModels$ff_file%d',i))) ##ok<NASGU> 
        strffModelFilePath = returnArgs[[0]]
        strffModelFilename  = returnArgs[[1]]
        ##         if (~strcmp(strffModelFilePath,''))
        ##             cd(strffModelFilePath);
        ##         end
        poped.db["model"]subffPointers[paste('ff_pointer',i,sep='')] = strffModelFilename
        i = i+1

	
      
    
        







    








			

	
	