"""
Author: Caiya Zhang, Yuchen Zheng
"""


import path
import random
import datetime
import numpy as np
import pandas as pd
from os import name
from project.sfg import sfg
from project.cell import cell
from project.size import size
from matpy.num import num
from matpy.matrix import Matrix
from project.zeros import zeros
from project.feval import feval
from project.pargen import pargen
from project.util import is_not_none
from project.getfulld import getfulld
from project.fileparts import fileparts
from project.poped_choose import poped_choose
from project.param_choose import param_choose
from project.param_set import param_set
from project.test_mat_size import test_mat_size
from project.create_design import create_design
from project.create_design_space import create_design_space
from project.find_largest_index import find_largest_index
from project.convert_variables import convert_variables
from project.get_all_params import get_all_params
from project.get_unfixed_params import get_unfixed_params


def reorder_vec(your_vec: Matrix, name_order):
    if your_vec.get_datanam() is not None:
        if all(your_vec.get_datanam()) in name_order:
            your_vec = your_vec[name_order[name_order in your_vec.get_datanam()]]

    return your_vec


def create_poped_database(pypkpdInput={},
                            ff_file=None,
                            ff_fun=None,
                            fg_file=None,
                            fg_fun=None,
                            fError_file=None,
                            fError_fun=None,
                            optsw=None,
                            # design,
                            # design_space,
                            iFIMCalculationType=None,
                            iApproximationMethod=None,
                            iFOCENumInd=None,
                            prior_fim=None,
                            strAutoCorrelationFile=None,
                            d_switch=None,
                            ofv_calc_type=None,
                            ds_index=None,
                            strEDPenaltyFile=None,
                            ofv_fun=None,
                            iEDCalculationType=None,
                            ED_samp_size=None,
                            bLHS=None,
                            strUserDistributionFile=None,
                            nbpop=None,
                            NumRanEff=None,
                            NumDocc=None,
                            NumOcc=None,
                            bpop=None,
                            d=None,
                            covd=None,
                            sigma=None,
                            docc=None,
                            covdocc=None,
                            notfixed_bpop=None,
                            notfixed_d=None,
                            notfixed_covd=None,
                            notfixed_docc=None,
                            notfixed_covdocc=None,
                            notfixed_sigma=None,
                            notfixed_covsigma=None,
                            bUseRandomSearch=None,
                            bUseStochasticGradient=None,
                            bUseLineSearch=None,
                            bUseExchangeAlgorithm=None,
                            bUseBFGSMinimizer=None,
                            EACriteria=None,
                            strRunFile=None,
                            pypkpd_version=None,
                            modtit=None,
                            output_file=None,
                            output_function_file=None,
                            strIterationFileName=None,
                            user_data=None,
                            ourzero=None,
                            dSeed=None,
                            line_opta=None,
                            line_optx=None,
                            bShowGraphs=None,
                            use_logfile=None,
                            m1_switch=None,
                            m2_switch=None,
                            hle_switch=None,
                            gradff_switch=None,
                            gradfg_switch=None,
                            grad_all_switch=None,
                            rsit_output=None,
                            sgit_output=None,
                            hm1=None,
                            hlf=None,
                            hlg=None,
                            hm2=None,
                            hgd=None,
                            hle=None,
                            AbsTol=None,
                            RelTol=None,
                            iDiffSolverMethod=None,
                            bUseMemorySolver=None,
                            rsit=None,
                            sgit=None,
                            intrsit=None,
                            intsgit=None,
                            maxrsnullit=None,
                            convergence_eps=None,
                            rslxt=None,
                            rsla=None,
                            cfaxt=None,
                            cfaa=None,
                            bGreedyGroupOpt=None,
                            EAStepSize=None,
                            EANumPoints=None,
                            EAConvergenceCriteria=None,
                            bEANoReplicates=None,
                            BFGSConvergenceCriteriaMinStep=None,
                            BFGSProjectedGradientTol=None,
                            BFGSTolerancef=None,
                            BFGSToleranceg=None,
                            BFGSTolerancex=None,
                            ED_diff_it=None,
                            ED_diff_percent=None,
                            line_search_it=None,
                            Doptim_iter=None,
                            iCompileOption=None,
                            iUseParallelMethod=None,
                            MCC_Dep=None,
                            strExecuteName=None,
                            iNumProcesses=None,
                            iNumChunkDesignEvals=None,
                            Mat_Out_Pre=None,
                            strExtraRunOptions=None,
                            dPollResultTime=None,
                            strFunctionInputName=None,
                            bParallelRS=None,
                            bParallelSG=None,
                            bParallelMFEA=None,
                            bParallelLS=None):

    # Filname and path of the model file
    ff_fun = param_choose(ff_fun, pypkpdInput, None, "model", "ff_pointer")
    fg_fun = param_choose(fg_fun, pypkpdInput, None, "model", "fg_pointer")
    fError_fun = param_choose(fError_fun, pypkpdInput, None, "model", "ferror_pointer")
    optsw = param_choose(optsw, pypkpdInput, Matrix(np.array([0, 0, 0, 0, 0])), "settings", "optsw")
    
    # Initial Design
    xt = param_choose(xt, pypkpdInput, Exception("'xt' needs to be defined"), "design", "xt")
    m = param_choose(m, pypkpdInput, None, "design", "m")
    x = param_choose(x, pypkpdInput, None, "design", "x")
    nx = param_choose(nx, pypkpdInput, None, "design", "nx")
    a = param_choose(a, pypkpdInput, None, "design", "a")
    groupsize = param_choose(groupsize, pypkpdInput, Exception("'groupsize' needs to be defined"), "design", "groupsize")
    ni = param_choose(ni, pypkpdInput, None, "design", "ni")
    model_switch = param_choose(model_switch, pypkpdInput, None, "design", "model_switch")
    
    # Design space
    maxni = param_choose(maxni, pypkpdInput, None, "design_space", "maxni")
    minni = param_choose(minni, pypkpdInput, None, "design_space", "minni")
    maxtotni = param_choose(maxtotni, pypkpdInput, None, "design_space", "maxtotni")
    mintotni = param_choose(mintotni, pypkpdInput, None, "design_space", "mintotni")
    maxgroupsize = param_choose(maxgroupsize, pypkpdInput, None, "design_space", "maxgroupsize")
    mingroupsize = param_choose(mingroupsize, pypkpdInput, None, "design_space", "mingroupsize")
    maxtotgroupsize = param_choose(maxtotgroupsize, pypkpdInput, None, "design_space", "maxtotgroupsize")
    mintotgroupsize = param_choose(mintotgroupsize, pypkpdInput, None, "design_space", "mintotgroupsize")
    maxxt = param_choose(maxxt, pypkpdInput, None, "design_space", "maxxt")
    minxt = param_choose(minxt, pypkpdInput, None, "design_space", "minxt")
    discrete_xt = param_choose(discrete_xt, pypkpdInput, None, "design_space", "xt_space")
    discrete_x = param_choose(discrete_x, pypkpdInput, None, "design_space", "discrete_x")
    maxa = param_choose(maxa, pypkpdInput, None, "design_space", "maxa")
    mina = param_choose(mina, pypkpdInput, None, "design_space", "mina")
    discrete_a = param_choose(discrete_a, pypkpdInput, None, "design_space", "a_space")
    bUseGrouped_xt = param_choose(bUseGrouped_xt, pypkpdInput, False, "design_space", "bUseGrouped_xt")
    G_xt = param_choose(G_xt, pypkpdInput, None, "design_space", "G_xt")
    bUseGrouped_a = param_choose(bUseGrouped_a, pypkpdInput, False, "design_space", "bUseGrouped_a")
    G_a = param_choose(G_a, pypkpdInput, None, "design_space", "G_a")
    bUseGrouped_x = param_choose(bUseGrouped_x, pypkpdInput, False, "design_space", "bUseGrouped_x")
    G_x = param_choose(G_x, pypkpdInput, None, "design_space", "G_x")

    # FIM calculation
    iFIMCalculationType = param_choose(iFIMCalculationType, pypkpdInput, 1, "settings", "iFIMCalculationType")
    iApproximationMethod = param_choose(iApproximationMethod, pypkpdInput, 0, "settings", "iApproximationMethod")
    iFOCENumInd = param_choose(iFOCENumInd, pypkpdInput, 1000, "settings", "iFOCENumInd")
    prior_fim = param_choose(prior_fim, pypkpdInput, Matrix(np.array([0, 0, 1])), "settings", "prior_fim")
    strAutoCorrelationFile = param_choose(strAutoCorrelationFile, pypkpdInput, "", "model", "auto_pointer")

    # Criterion specification
    d_switch = param_choose(d_switch, pypkpdInput, 1, "settings", "d_switch")
    ofv_calc_type = param_choose(ofv_calc_type, pypkpdInput, 4, "settings", "ofv_calc_type")
    ds_index = param_set(ds_index, pypkpdInput, "parameters", "ds_index")
    strEDPenaltyFile = param_choose(strEDPenaltyFile, pypkpdInput, "", "settings", "strEDPenaltyFile")
    ofv_fun = param_choose(ofv_fun, pypkpdInput, None, "settings", "ofv_fun")
    
    # E-family Criterion options
    iEDCalculationType = param_choose(iEDCalculationType, pypkpdInput, 0, "settings", "iEDCalculationType")
    ED_samp_size = param_choose(ED_samp_size, pypkpdInput, 45, "settings", "ED_samp_size")
    bLHS = param_choose(bLHS, pypkpdInput, 1, "settings", "bLHS")
    strUserDistributionFile = param_choose(strUserDistributionFile, pypkpdInput, "", "model", "user_distribution_pointer")

    # Model parameters
    nbpop = param_set(ds_index, pypkpdInput, "parameters", "nbpop")
    NumRanEff = param_set(NumRanEff, pypkpdInput, "parameters", "NumRanEff")
    NumDocc = param_set(NumDocc, pypkpdInput, "parameters", "NumDocc")
    NumOcc = param_set(NumOcc, pypkpdInput, "parameters", "NumOcc")
    bpop = param_choose(bpop, pypkpdInput, Exception("bpop must be defined"), "parameters", "bpop")
    d = param_choose(d, pypkpdInput, None, "parameters", "d")
    covd = param_set(covd, pypkpdInput, "parameters", "covd")
    sigma = param_set(sigma, pypkpdInput, "parameters", "sigma")
    docc = param_choose(docc, pypkpdInput, Matrix(np.array([0, 0, 3])), "parameters", "docc")
    covdocc = param_choose(covdocc, pypkpdInput, np.zeros([1, docc.get_shape()[1]*(docc.get_shape()[1]-1)/2]), "parameters", "covdocc")

    # Model parameters fixed or not
    notfixed_bpop = param_set(notfixed_bpop, pypkpdInput, "parameters", "notfixed_bpop")
    notfixed_d = param_set(notfixed_d, pypkpdInput, "parameters", "notfixed_d")
    notfixed_covd = param_set(notfixed_covd, pypkpdInput, "parameters", "notfixed_covd")
    notfixed_docc = param_set(notfixed_docc, pypkpdInput, "parameters", "notfixed_docc")
    notfixed_covdocc = param_choose(notfixed_covdocc, pypkpdInput, Matrix(np.zeros([1, covdocc.get_size()])), "parameters", "notfixed_covdocc")
    notfixed_sigma = param_choose(notfixed_sigma, pypkpdInput, Matrix([1]*sigma.get_shape()[1]), "parameters", "notfixed_sigma")
    notfixed_covsigma = param_choose(notfixed_covsigma, pypkpdInput, Matrix(np.zeros([1, notfixed_sigma.get_size()*(notfixed_sigma.get_size()-1)/2])), "parameters", "notfixed_covsigma")

    # Optimization algorithm choices
    bUseRandomSearch = param_choose(bUseRandomSearch, pypkpdInput, True, "settings", "bUseRandomSearch")
    bUseStochasticGradient = param_choose(bUseStochasticGradient, pypkpdInput, True, "settings", "bUseStochasticGradient")
    bUseLineSearch = param_choose(bUseLineSearch, pypkpdInput, True, "settings", "bUseLineSearch")
    bUseExchangeAlgorithm = param_choose(bUseExchangeAlgorithm, pypkpdInput, False, "settings", "bUseExchangeAlgorithm")
    bUseBFGSMinimizer = param_choose(bUseBFGSMinimizer, pypkpdInput, False, "settings", "bUseBFGSMinimizer")
    EACriteria = param_choose(EACriteria, pypkpdInput, 1, "settings", "EACriteria")
    strRunFile = param_choose(strRunFile, pypkpdInput, "", "settings", "run_file_pointer")

    # Labeling and file names
    pypkpd_version = param_choose(pypkpd_version, pypkpdInput, pypkpd.__version__, "settings", "pypkpd_version")
    modtit = param_choose(modtit, pypkpdInput, "pypkpd model", "settings", "modtit")
    output_file = param_choose(output_file, pypkpdInput, "pypkpd_output_summary", "settings", "output_file")
    output_function_file = param_choose(output_function_file, pypkpdInput, "pypkpd_output", "settings", "output_function_file")
    strIterationFileName = param_choose(strIterationFileName, pypkpdInput, "pypkpd_current.py", "settings", "strIterationFileName")
    
    # Misc options

    # Parallel options for pypkpd


#####-------------- main part --------------#####
    pypkpd_db = {}

    pypkpd_db["settings"] = {}
    pypkpd_db["settings"]["pypkpd_version"] = pypkpd_version

    # if BFGSConvergenceCriteriaMinStep is None:
    #     BFGSConvergenceCriteriaMinStep = param_choose(pypkpdInput, 1e-08, 0, "settings", "BFGSConvergenceCriteriaMinStep")

    # if MCC_Dep is None:
    #     MCC_Dep = param_choose(pypkpdInput, "", 0, "settings", "parallel", "strAdditionalMCCCompilerDependencies")

    pypkpd_db["model"] = {}
    pypkpd_db["model"]["user_distribution_pointer"] = ""

    if strUserDistributionFile is not None:
        if str(strUserDistributionFile) != "":
            pypkpd_db["model"]["user_distribution_pointer"] = strUserDistributionFile
    else:
        exec(open(strUserDistributionFile).read())
        returnArgs = fileparts(strUserDistributionFile)
        strUserDistFilePath = list(returnArgs.values())[0]
        strUserDistFilename = list(returnArgs.values())[1]
        pypkpd_db["model"]["user_distribution_pointer"] = strUserDistFilename

    # should be removed
    if (np.array(x.get_shape()) == 0).any():
        x = None
        G_x = None
        discrete_x = None

    # should be removed
    if (np.array(a.get_shape()) == 0).any():
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

    # design_space["design"] = design
    # design_space["design_space"] = design_space

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

    # should be removed
    if ~is_not_none(design, "x"):
        design["x"] = Matrix(np.zeros([design["m"].get_value(), 0]))
        design_space["G_x"] = design["x"]
        design_space["bUseGrouped_x"] = False
        design_space["discrete_x"] = cell([design["m"].get_value(), 0])

    # should be removed
    if ~is_not_none(design, "a"):
        design["a"] = Matrix(np.zeros([design["m"].get_value(), 0]))
        design_space["G_a"] = design["a"]
        design_space["bUseGrouped_a"] = False
        design_space["mina"] = design["a"]
        design_space["maxa"] = design["a"]

    pypkpd_db["design"] = design
    pypkpd_db["design_space"] = design_space

    pypkpd_db["settings"] = {}
    pypkpd_db["settings"]["bLHS"] = bLHS

# pypkpd_db$discrete_x = design_space$discrete_x # should be removed only in design_space

# pypkpd_db$maxni=max(design_space$maxni) # should be only in design_space and called maxmaxni if needed
# pypkpd_db$minni=min(design_space$minni) # should be only in design_space and called minminni if needed

# pypkpd_db$bUseGrouped_xt = design_space$bUseGrouped_xt # should be only in design_space
# pypkpd_db$bUseGrouped_a  = design_space$bUseGrouped_a # should be only in design_space
# pypkpd_db$bUseGrouped_x  = design_space$bUseGrouped_x # should be only in design_space

    pypkpd_db["settings"]["d_switch"] = d_switch
    pypkpd_db["settings"]["iApproximationMethod"] = iApproximationMethod

    pypkpd_db["settings"]["iFOCENumInd"] = iFOCENumInd
    pypkpd_db["settings"]["bUseRandomSearch"] = bUseRandomSearch
    pypkpd_db["settings"]["bUseStochasticGradient"] = bUseStochasticGradient
    pypkpd_db["settings"]["bUseLineSearch"] = bUseLineSearch
    pypkpd_db["settings"]["bUseExchangeAlgorithm"] = bUseExchangeAlgorithm
    pypkpd_db["settings"]["bUseBFGSMinimizer"] = bUseBFGSMinimizer

    pypkpd_db["settings"]["iEDCalculationType"] = iEDCalculationType

    pypkpd_db["settings"]["BFGSConvergenceCriteriaMinStep"] = BFGSConvergenceCriteriaMinStep
    pypkpd_db["settings"]["BFGSProjectedGradientTol"] = BFGSProjectedGradientTol
    pypkpd_db["settings"]["BFGSTolerancef"] = BFGSTolerancef
    pypkpd_db["settings"]["BFGSToleranceg"] = BFGSToleranceg
    pypkpd_db["settings"]["BFGSTolerancex"] = BFGSTolerancex

    pypkpd_db["parameters"] = {}
    pypkpd_db["parameters"]["covdocc"] = covdocc
    pypkpd_db["parameters"]["notfixed_covdocc"] = notfixed_covdocc
    pypkpd_db["parameters"]["notfixed_covsigma"] = notfixed_covsigma

    pypkpd_db["settings"]["parallel"] = {}
    pypkpd_db["settings"]["parallel"]["iCompileOption"] = iCompileOption
    pypkpd_db["settings"]["parallel"]["strAdditionalMCCCompilerDependencies"] = MCC_Dep
    pypkpd_db["settings"]["parallel"]["iUseParallelMethod"] = iUseParallelMethod
    pypkpd_db["settings"]["parallel"]["strExecuteName"] = strExecuteName
    pypkpd_db["settings"]["parallel"]["iNumProcesses"] = iNumProcesses
    pypkpd_db["settings"]["parallel"]["iNumChunkDesignEvals"] = iNumChunkDesignEvals
    #pypkpd_db["settings"]["parallel"]["strMatFileInputPrefix"] = strMatFileInputPrefix
    pypkpd_db["settings"]["parallel"]["strMatFileOutputPrefix"] = Mat_Out_Pre
    pypkpd_db["settings"]["parallel"]["strExtraRunOptions"] = strExtraRunOptions
    pypkpd_db["settings"]["parallel"]["dPollResultTime"] = dPollResultTime
    pypkpd_db["settings"]["parallel"]["strFunctionInputName"] = strFunctionInputName
    pypkpd_db["settings"]["parallel"]["bParallelRS"] = bParallelRS
    pypkpd_db["settings"]["parallel"]["bParallelSG"] = bParallelSG
    pypkpd_db["settings"]["parallel"]["bParallelLS"] = bParallelLS
    pypkpd_db["settings"]["parallel"]["bParallelMFEA"] = bParallelMFEA

    pypkpd_db["settings"]["hm1"] = hm1
    pypkpd_db["settings"]["hlf"] = hlf
    pypkpd_db["settings"]["hlg"] = hlg
    pypkpd_db["settings"]["hm2"] = hm2
    pypkpd_db["settings"]["hgd"] = hgd
    pypkpd_db["settings"]["hle"] = hle

    pypkpd_db["settings"]["AbsTol"] = AbsTol
    pypkpd_db["settings"]["RelTol"] = RelTol
    pypkpd_db["settings"]["iDiffSolverMethod"] = iDiffSolverMethod
    # Temp thing for memory solvers
    pypkpd_db["settings"]["bUseMemorySolver"] = bUseMemorySolver
    pypkpd_db["settings"]["solved_solutions"] = cell(0, 0)
    pypkpd_db["settings"]["maxtime"] = str(max(max(pypkpd_db["design_space"]["maxxt"].get_all_data()))) + str(pypkpd_db["settings"]["hgd"])

    pypkpd_db["settings"]["iFIMCalculationType"] = iFIMCalculationType
    pypkpd_db["settings"]["rsit"] = rsit
    pypkpd_db["settings"]["sgit"] = sgit
    pypkpd_db["settings"]["intrsit"] = intrsit
    pypkpd_db["settings"]["intsgit"] = intsgit
    pypkpd_db["settings"]["maxrsnullit"] = maxrsnullit
    pypkpd_db["settings"]["convergence_eps"] = convergence_eps
    pypkpd_db["settings"]["rslxt"] = rslxt
    pypkpd_db["settings"]["rsla"] = rsla
    pypkpd_db["settings"]["cfaxt"] = cfaxt
    pypkpd_db["settings"]["cfaa"] = cfaa

    pypkpd_db["settings"]["EACriteria"] = EACriteria
    pypkpd_db["settings"]["EAStepSize"] = EAStepSize
    pypkpd_db["settings"]["EANumPoints"] = EANumPoints
    pypkpd_db["settings"]["EAConvergenceCriteria"] = EAConvergenceCriteria

    pypkpd_db["settings"]["ED_samp_size"] = ED_samp_size
    pypkpd_db["settings"]["ED_diff_it"] = ED_diff_it
    pypkpd_db["settings"]["ED_diff_percent"] = ED_diff_percent
    pypkpd_db["settings"]["ls_step_size"] = line_search_it

    pypkpd_db["settings"]["ofv_calc_type"] = ofv_calc_type

    pypkpd_db["settings"]["iNumSearchIterationsIfNotLineSearch"] = Doptim_iter

    pypkpd_db["settings"]["ourzero"] = ourzero
    pypkpd_db["settings"]["rsit_output"] = rsit_output
    pypkpd_db["settings"]["sgit_output"] = sgit_output

    if callable(fg_fun):
        pypkpd_db["model"]["fg_pointer"] = fg_fun
    elif type(fg_fun) is str:
        if fg_fun is not None:
            pypkpd_db["model"]["fg_pointer"] = fg_fun
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
            pypkpd_db["model"]["fg_pointer"] = strfgModelFilename
        else:
            pypkpd_db["model"]["fg_pointer"] = fg_file

    pypkpd_db["settings"]["ed_penalty_pointer"] = zeros(1, 0)
    if str(strEDPenaltyFile) != "":
        if strEDPenaltyFile is not None:
            pypkpd_db["settings"]["ed_penalty_pointer"] = strEDPenaltyFile
        else:
            exec(open(pypkpdInput["strEDPenaltyFile"]).read())
            returnArgs = fileparts(pypkpdInput["strEDPenaltyFile"])
            strEDPenaltyFilePath = list(returnArgs.values())[0]
            strEDPenaltyFilename = list(returnArgs.values())[1]
    # if (~strcmp(strEDPenaltyFilePath,''))
    # cd(strEDPenaltyFilePath);
    # end
            pypkpd_db["settings"]["ed_penalty_pointer"] = strEDPenaltyFilename

    # if(is.null(ofv_fun) or is.function(ofv_fun)){
#   pypkpd_db["settings"]ofv_fun = ofv_fun
# } else {
#   stop("ofv_fun must be a function or None")
# }

    if ofv_fun is None or callable(ofv_fun):
        pypkpd_db["settings"]["ofv_fun"] = ofv_fun
    else:
        # source explicit file
        # here I assume that function in file has same name as filename minus .txt and pathnames
        if ofv_fun is not None:
            exec(open(str(ofv_fun)).read())
            pypkpd_db["settings"]["ofv_fun"] = eval('text=fileparts(ofv_fun)["filename"]')
        else:
            raise Exception("ofv_fun is not a function or None, and no file with that name was found")

    # if(is.function(ofv_fun)){
#   pypkpd_db["settings"]ofv_fun = ofv_fun
# } else if(exists(ofv_fun)){
#   pypkpd_db["settings"]ofv_fun = ofv_fun
# } else {
#   source(ofv_fun)
#   returnArgs =  fileparts(ofv_fun)
#   strffModelFilePath = returnArgs[1]
#   strffModelFilename  = returnArgs[1]
#   ## if (~strcmp(strffModelFilePath,''))
#   ##    cd(strffModelFilePath);
#   ## end
#   pypkpd_db["settings"]ofv_fun = strffModelFilename
# }

    pypkpd_db["model"]["auto_pointer"] = zeros(1, 0)
    #pypkpd_db["model"]["auto_pointer"] = ''
    if strAutoCorrelationFile != "":
        if str(strAutoCorrelationFile) != "":
            if strAutoCorrelationFile is not None:
                pypkpd_db["model"]["auto_pointer"] = strAutoCorrelationFile
            else:
                exec(open(pypkpdInput["strAutoCorrelationFile"]).read())
                returnArgs = fileparts(pypkpdInput["strAutoCorrelationFile"])
                strAutoCorrelationFilePath = list(returnArgs.values())[0]
                strAutoCorrelationFilename = list(returnArgs.values())[1]
            # if (~strcmp(strAutoCorrelationFilePath,''))
            # cd(strAutoCorrelationFilePath);
            # end
                pypkpd_db["model"]["auto_pointer"] = strAutoCorrelationFilename

    if callable(ff_fun):
        pypkpd_db["model"]["ff_pointer"] = ff_fun
    elif type(ff_fun) is str:
        if ff_fun is not None:
            pypkpd_db["model"]["ff_pointer"] = ff_fun
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
            pypkpd_db["model"]["ff_pointer"] = strffModelFilename
        else:
            pypkpd_db["model"]["ff_pointer"] = ff_file

    # Check if there is any sub models defined
    for key in pypkpdInput:
        if "SubModels" in pypkpdInput.keys():
            i = 1
            for elem_key in pypkpdInput["SubModels"]:
                if elem_key == ("ff_file%d" % i):
                    # ok<np.nanSGU>
                    exec(eval('pypkpdInput["SubModels"]["ff_file"]%d' % i))
                    returnArgs = fileparts(eval('pypkpdInput["SubModels"]["ff_file"]%d' % i))  # ok<np.nanSGU>
                    strffModelFilePath = list(returnArgs.values())[0]
                    strffModelFilename = list(returnArgs.values())[1]
                    # if (~strcmp(strffModelFilePath,''))
                    # cd(strffModelFilePath);
                    # end
                    pypkpd_db["model"]["subffPointers"]["".join(["ff_pointer", str(i)])] = strffModelFilename
                    i = i + 1

    if callable(fError_fun):
        pypkpd_db["model"]["ferror_pointer"] = fError_fun
    elif type(fError_fun) is str:
        if fError_fun is not None:
            pypkpd_db["model"]["ferror_pointer"] = fError_fun
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
            pypkpd_db["model"]["ferror_pointer"] = strErrorModelFilename
        else:
            pypkpd_db["model"]["ferror_pointer"] = fError_file

    # %Set the model file string path
# pypkpd_db$model_file = ff_file
##   model_file = eval('functions(pypkpd_db.ff_pointer)');
# if (~strcmp(model_file.file,''))
##       pypkpd_db.model_file = eval('char(model_file.file)');
# else
##       pypkpd_db.model_file = eval('char(model_file.function)');
# end

# ==================================
# Initialize the randomization
# ==================================
    if dSeed is not None:
        if dSeed == -1:
            pypkpd_db["settings"]["dSeed"] = datetime.datetime.now()
        else:
            pypkpd_db["settings"]["dSeed"] = dSeed
        random.seed(pypkpd_db["settings"]["dSeed"])

    pypkpd_db["parameters"]["nbpop"] = poped_choose(nbpop, find_largest_index(pypkpd_db["model"]["fg_pointer"], "bpop"), 0)
    pypkpd_db["parameters"]["NumRanEff"] = poped_choose(NumRanEff, find_largest_index(pypkpd_db["model"]["fg_pointer"], "b"), 0)
    pypkpd_db["parameters"]["NumDocc"] = poped_choose(NumDocc, find_largest_index(pypkpd_db["model"]["fg_pointer"], "bocc", mat=True, mat_row=True), 0)
    pypkpd_db["parameters"]["NumOcc"] = poped_choose(NumOcc, find_largest_index(pypkpd_db["model"]["fg_pointer"], "bocc", mat=True, mat_row=False), 0)
    
    pypkpd_db["parameters"]["ng"] = eval(str(pypkpd_db["model"]["fg_pointer"]) + "(num(0), num(0), num(0), num(0)," + str(zeros(pypkpd_db["parameters"]["NumDocc"], pypkpd_db["parameters"]["NumOcc"])) + ")").get_size()
    
    docc_arr = np.array([1])
    d_arr = np.array([1])
    bpop_arr = np.array([1])
    pypkpd_db["parameters"]["notfixed_docc"] = poped_choose(notfixed_docc, 
        Matrix(np.ones([1, param_choose(pypkpd_db, 0, 0, None, None, "parameters", "NumDocc")])), 0)
    pypkpd_db["parameters"]["notfixed_d"] = poped_choose(notfixed_d, 
        Matrix(np.ones([1, param_choose(pypkpd_db, 0, 0, None, None, "parameters", "NumRanEff")])), 0)
    pypkpd_db["parameters"]["notfixed_bpop"] = poped_choose(notfixed_bpop, 
        Matrix(np.ones([1, param_choose(pypkpd_db, 0, 0, None, None, "parameters", "nbpop")])), 0)


# reorder named values
    fg_names = eval(pypkpd_db["model"]["fg_pointer"] + "(num(1), num(1), num(1), num(1)," + str(ones(param_choose(pypkpd_db, 0, 0, None, None, "parameters", "NumDocc"), param_choose(pypkpd_db, 0, 0, None, None, "parameters", "NumDocc"))) + ")").get_datanam()
    
    pypkpd_db["parameters"]["notfixed_bpop"] = reorder_vec(pypkpd_db["parameters"]["notfixed_bpop"], fg_names)

    pypkpd_db["parameters"]["notfixed_d"] = reorder_vec(pypkpd_db["parameters"]["notfixed_d"], fg_names)

    # we have just the parameter values not the uncertainty
    if size(d)[0] == 1 and size(d)[1] == pypkpd_db["parameters"]["NumRanEff"]:
        d_descr = zeros(pypkpd_db["parameters"]["NumRanEff"], 3)

        # reorder named values
        d = reorder_vec(d, fg_names)
        
        #set_data
        #d_descr.set_data()[:,1] = d
        for i in range(0, d.get_size()):
            d_descr.set_one_data(d[i], [i,1], [i,1])
        d_descr.set_multiple_data(0, [0,d_descr.get_shape[0]], [0,1])  # point values
        d_descr.set_multiple_data(0, [0,d_descr.get_shape[0]], [2,3])  # variance
        d_descr.set_colnam(d.get_datanam())
        # d_descr.columns.values = d.keys()
        d = d_descr

    # we have just the parameter values not the uncertainty
    if size(bpop)[0] == 1 and size(bpop)[1] == pypkpd_db["parameters"]["nbpop"]:
        bpop_descr = zeros(pypkpd_db["parameters"]["nbpop"], 3)

        # reorder named values
        bpop = reorder_vec(bpop, fg_names)

        
        #set_data
        #bpop_descr.set_data()[:,1] = bpop
        for i in range(0, bpop.get_size()):
            bpop_descr.set_one_data(bpop[i], [i,1], [i,1])
        bpop_descr.set_multiple_data(0, [0,bpop_descr.get_shape[0]], [0,1])  # point values
        bpop_descr.set_multiple_data(0, [0,bpop_descr.get_shape[0]], [2,3])  # variance
        bpop_descr.set_colnam(bpop.get_datanam())
        # bpop_descr.columns.values = bpop.keys()
        bpop = bpop_descr

    # we have just the diagonal parameter values
    if size(sigma)[0] == 1 and (type(sigma) is not np.ndarray or type(sigma) is not Matrix):
        sigma_tmp = Matrix(np.diag(sigma, size(sigma)[1], size(sigma)[1]))
        sigma_tmp.set_colnam(sigma.keys())
        # sigma_tmp.columns.values = sigma.keys()
        sigma = sigma_tmp

    covd = poped_choose(covd, zeros(
        1, pypkpd_db["parameters"]["NumRanEff"])*(pypkpd_db["parameters"]["NumRanEff"]-1)/2, 0)
    pypkpd_db["parameters"]["covd"] = covd

    tmp = ones(1, len(covd))
    if tmp is not None:
        for i in range(0, len(covd)):
            if covd[i] == 0:
                tmp[i] = 0
    #tmp[covd==0] = 0
    pypkpd_db["parameters"]["notfixed_covd"] = poped_choose(notfixed_covd, tmp, 0)

    # ==================================
    # Sample the individual eta's for FOCE and FOCEI
    # ==================================
    if pypkpd_db["settings"]["iApproximationMethod"] != 0 and pypkpd_db["settings"]["iApproximationMethod"] != 3:

        iMaxCorrIndNeeded = 100
        bzeros = zeros(pypkpd_db["parameters"]["NumRanEff"], 1)
        bones = Matrix(np.ones(int(pypkpd_db["parameters"]["NumRanEff"])), (int(pypkpd_db["parameters"]["NumRanEff"]), 1))
        bocczeros = zeros(pypkpd_db["parameters"]["NumDocc"], 1)
        boccones = Matrix(np.ones(int(pypkpd_db["parameters"]["NumDocc"])), (int(pypkpd_db["parameters"]["NumDocc"]), 1))

        pypkpd_db["parameters"]["b_global"] = zeros(pypkpd_db["parameters"]["NumRanEff"], max(
            pypkpd_db["settings"]["iFOCENumInd"], iMaxCorrIndNeeded))

        fulld = getfulld(Matrix(d.get_all_data()[:,1]), pypkpd_db["parameters"]["covd"])
        fulldocc = getfulld(Matrix(docc.get_all_data()[:,1]), pypkpd_db["parameters"]["covdocc"])

        pypkpd_db["parameters"]["bocc_global"] = cell(
            pypkpd_db["settings"]["iFOCENumInd"], 1)

        if pypkpd_db["settings"]["d_switch"] is True:
            pypkpd_db["parameters"]["b_global"] = Matrix(np.transpose(np.random.multivariate_normal(
                max(pypkpd_db["settings"]["iFOCENumInd"], iMaxCorrIndNeeded), sigma=fulld)))
            for i in range(0, pypkpd_db["settings"]["iFOCENumInd"]):
                pypkpd_db["parameters"]["bocc_global"][i] = zeros(
                    size(docc)[0], pypkpd_db["parameters"]["NumOcc"])
                if pypkpd_db["parameters"]["NumOcc"] != 0:
                    pypkpd_db["parameters"]["bocc_global"][i] = Matrix(np.transpose(np.random.multivariate_normal(
                        pypkpd_db["parameters"]["NumOcc"], sigma=fulldocc)))

        else:
            d_dist = pargen(d, pypkpd_db["model"]["user_distribution_pointer"], max(
                pypkpd_db["settings"]["iFOCENumInd"], iMaxCorrIndNeeded), pypkpd_db["settings"]["bLHS"], zeros(1, 0), pypkpd_db)
            docc_dist = pargen(docc, pypkpd_db["model"]["user_distribution_pointer"], pypkpd_db["settings"]
                               ["iFOCENumInd"], pypkpd_db["settings"]["bLHS"], zeros(1, 0), pypkpd_db)

            ## set one data with type of Matrix
            if len(d_dist) != 0:
                for i in range(0, max(pypkpd_db["settings"]["iFOCENumInd"], iMaxCorrIndNeeded)):
                    for j in range(0, pypkpd_db["parameters"]["b_global"].get_shape[0]):
                        pypkpd_db["parameters"]["b_global"].set_one_data(Matrix(np.transpose(
                            np.random.multivariate_normal(
                                1, sigma=getfulld(d_dist[i,:], pypkpd_db["parameters"]["covd"])))), None, [j,i])                        

            if len(docc_dist) != 0:
                for i in range(0, pypkpd_db["settings"]["iFOCENumInd"]):
                    pypkpd_db["parameters"]["bocc_global"].set_one__data(Matrix(np.transpose(
                        np.random.multivariate_normal(
                        pypkpd_db["parameters"]["NumOcc"], sigma=getfulld(docc_dist[i,:], pypkpd_db["parameters"]["covdocc"])))), None, [0,i])
    else:
        pypkpd_db["parameters"]["bocc_global"] = zeros(pypkpd_db["parameters"]["NumRanEff"], 1)
        pypkpd_db["parameters"]["bocc_global"] = cell(1, 1)
        pypkpd_db["parameters"]["bocc_global"].set_one_data(zeros(size(docc)[0], pypkpd_db["parameters"]["NumOcc"]), None, [0,1])
        
        pypkpd_db["settings"]["iFOCENumInd"] = 1

    pypkpd_db["settings"]["modtit"] = modtit
    pypkpd_db["settings"]["exptit"] = ('%s_exp["mat"]', modtit)  # experiment settings title
    pypkpd_db["settings"]["opttit"] = ('%s_opt["mat"]', modtit)  # optimization settings title
    pypkpd_db["settings"]["bShowGraphs"] = bShowGraphs

    pypkpd_db["settings"]["use_logfile"] = use_logfile
    pypkpd_db["settings"]["output_file"] = output_file
    pypkpd_db["settings"]["output_function_file"] = output_function_file

    pypkpd_db["settings"]["optsw"] = optsw

    line_opta = poped_choose(line_opta, np.ones([1, size(pypkpd_db["design"]["a"])[1]]), 0)
    if test_mat_size(np.array([1, size(pypkpd_db["design"]["a"])[0]]), np.array(line_opta), "line_opta"):
        pypkpd_db["settings"]["line_opta"] = line_opta

    line_optx = poped_choose(line_optx, np.ones([1, size(pypkpd_db["design"]["x"])[1]]), 0)
    if test_mat_size(np.array([1, size(pypkpd_db["design"]["x"])[0]]), np.array(line_opta), "line_optx"):
        pypkpd_db["settings"]["line_optx"] = line_optx


# pypkpd_db$design = pypkpdInput$design
    pypkpd_db["parameters"]["bpop"] = bpop
    pypkpd_db["parameters"]["d"] = d
    pypkpd_db["parameters"]["covd"] = covd
    pypkpd_db["parameters"]["sigma"] = sigma
    pypkpd_db["parameters"]["docc"] = docc
    pypkpd_db["parameters"]["covdocc"] = covdocc

# pypkpd_db["design"]G = design_space$G_xt ## should only be in design_space
# pypkpd_db["design"]Ga = design_space$G_a ## should only be in design_space
# pypkpd_db["design"]Gx = design_space$G_x ## should only be in design_space

# pypkpd_db["design"]maxgroupsize = design_space$maxgroupsize ## should only be in design_space
# pypkpd_db["design"]mingroupsize = design_space$mingroupsize ## should only be in design_space
# pypkpd_db["design"]maxtotgroupsize = design_space$maxtotgroupsize ## should only be in design_space
# pypkpd_db["design"]mintotgroupsize = design_space$mintotgroupsize ## should only be in design_space
# pypkpd_db["design"]maxxt = design_space$maxxt ## should only be in design_space
# pypkpd_db["design"]minxt = design_space$minxt ## should only be in design_space
# pypkpd_db["design"]discrete_x = design_space$discrete_x ## should only be in design_space
# pypkpd_db["design"]maxa = design_space$maxa ## should only be in design_space
# pypkpd_db["design"]mina = design_space$mina ## should only be in design_space

    pypkpd_db["settings"]["m1_switch"] = m1_switch
    pypkpd_db["settings"]["m2_switch"] = m2_switch
    pypkpd_db["settings"]["hle_switch"] = hle_switch
    pypkpd_db["settings"]["gradff_switch"] = gradff_switch
    pypkpd_db["settings"]["gradfg_switch"] = gradfg_switch
    pypkpd_db["settings"]["grad_all_switch"] = grad_all_switch

    pypkpd_db["settings"]["prior_fim"] = prior_fim

    pypkpd_db["parameters"]["notfixed_sigma"] = notfixed_sigma
    # pypkpd_db["parameters"]sigma = sigma
    # pypkpd_db["parameters"]docc = docc

    pypkpd_db["parameters"]["ds_index"] = ds_index

    # create ds_index if not already done
    if pypkpd_db["parameters"]["ds_index"] is not None:
        unfixed_params = get_unfixed_params(pypkpd_db)
        pypkpd_db["parameters"]["ds_index"] = Matrix(np.transpose(
            np.repeat(0, unfixed_params["all"])))
        pypkpd_db["parameters"]["ds_index"][(
            unfixed_params["bpop"].get_size()+1):pypkpd_db["parameters"]["ds_index"].get_size()] = 1
    else:
        if type(pypkpd_db["parameters"]["ds_index"]) is not Matrix:
            pypkpd_db["parameters"]["ds_index"] = Matrix(np.array(
                [pypkpd_db["parameters"]["ds_index"], 1, len(pypkpd_db["parameters"]["ds_index"])]))

    pypkpd_db["settings"]["strIterationFileName"] = strIterationFileName
    pypkpd_db["settings"]["user_data"] = user_data
    pypkpd_db["settings"]["bUseSecondOrder"] = False
    pypkpd_db["settings"]["bCalculateEBE"] = False
    pypkpd_db["settings"]["bGreedyGroupOpt"] = bGreedyGroupOpt
    pypkpd_db["settings"]["bEANoReplicates"] = bEANoReplicates

    if len(strRunFile) != 0:
        if strRunFile == " ":
            pypkpd_db["settings"]["run_file_pointer"] = zeros(1, 0)
        else:
            if strRunFile is not None:
                pypkpd_db["settings"]["run_file_pointer"] = strRunFile
            else:
                exec(open(pypkpdInput["strRunFile"]).read())
                returnArgs = fileparts(pypkpdInput["strRunFile"])
                strRunFilePath = list(returnArgs.values())[0]
                strRunFilename = list(returnArgs.values())[1]
                # if (~strcmp(strErrorModelFilePath,''))
                # cd(strErrorModelFilePath);
                # end
                pypkpd_db["settings"]["run_file_pointer"] = strRunFilename

#pypkpd_db = convert_pypkpdInput(pypkpdInput,...)

    pypkpd_db["settings"]["Engine"] = {"Type": 1, "Version": poped_version["version"]}

    pypkpd_db = convert_variables(pypkpd_db)  # need to transform here

    param_val = get_all_params(pypkpd_db)
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
	sigma_full = pypkpd_db["parameters"]["sigma"]
	
	pypkpd_db["parameters"]["param_pt_val"]["bpop"] = bpop_val
	pypkpd_db["parameters"]["param_pt_val"]["d"] = d_full
	pypkpd_db["parameters"]["param_pt_val"]["docc"] = docc_full
	pypkpd_db["parameters"]["param_pt_val"]["sigma"] = sigma_full
    
	"""
    #     downsize.list = downsizing_general_design(pypkpd_db)
    #     tmp.names = names(downsize.list)
    #     model_switch = c()
    #     ni = c()
    #     xt = c()
    #     x = c()
    #     a = c()
    #     eval(parse(text=paste(tmp.names,"=","downsize.list$",tmp.names,sep="")))
    #     pypkpd_db$downsized.design$model_switch = model_switch
    #     pypkpd_db$downsized.design$ni = ni
    #     pypkpd_db$downsized.design$xt = xt
    #     pypkpd_db$downsized.design$x = x
    #     pypkpd_db$downsized.design$a = a
    #     pypkpd_db$downsized.design$groupsize = pypkpd_db["design"]groupsize

    retargs = fileparts(pypkpd_db["settings"]["output_file"])
    pypkpd_db["settings"]["strOutputFilePath"] = list(retargs.values())[0]
    pypkpd_db["settings"]["strOutputFileName"] = list(retargs.values())[1]
    pypkpd_db["settings"]["strOutputFileExtension"] = list(retargs.values())[2]

    return pypkpd_db

def somestring(**kwargs):
    return ", ".join(f"{key}={value}" for key, value in kwargs.items())