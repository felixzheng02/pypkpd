"""
Author: Caiya Zhang, Yuchen Zheng
"""


import path
import random
import datetime
import numpy as np
from project.cell import cell
from project.size import size
from project.length import length
from project.names import names
from matpy.num import Num
from matpy.matrix import Matrix
from project.pargen import pargen
from project.util import default_if_none, get_dict_value
from project.getfulld import getfulld
from project.fileparts import fileparts
from project.param_choose import param_choose
from project.pypkpd_choose import pypkpd_choose
from project.test_mat_size import test_mat_size
from project.create_design import create_design
from project.create_design_space import create_design_space
from project.find_largest_index import find_largest_index
from project.convert_variables import convert_variables
from project.get_all_params import get_all_params
from project.get_unfixed_params import get_unfixed_params
from project.call_func import call_func


def reorder_vec(your_vec: Matrix, name_order):
    if your_vec.get_datanam() is not None:
        if all(your_vec.get_datanam() in name_order):
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
                            xt = None,
                            m = None,
                            x = None,
                            nx = None,
                            a = None,
                            groupsize = None,
                            ni = None,
                            model_switch = None,
                            maxni = None,
                            minni = None,
                            maxtotni = None,
                            mintotni = None,
                            maxgroupsize = None,
                            mingroupsize = None,
                            maxtotgroupsize = None,
                            mintotgroupsize = None,
                            maxxt = None,  
                            minxt = None,
                            discrete_xt = None,
                            discrete_x = None,
                            maxa = None,
                            mina = None,
                            discrete_a = None,
                            bUseGrouped_xt = None,
                            G_xt = None,
                            bUseGrouped_a = None,
                            G_a = None,
                            bUseGrouped_x = None,
                            G_x = None,
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
    ff_fun = param_choose(ff_fun, get_dict_value(pypkpdInput, "model", "ff_pointer"), None)
    fg_fun = param_choose(fg_fun, get_dict_value(pypkpdInput, "model", "fg_pointer"), None)
    fError_fun = param_choose(fError_fun, get_dict_value(pypkpdInput, "model", "ferror_pointer"), None)
    optsw = param_choose(optsw, get_dict_value(pypkpdInput, "settings", "optsw"), np.zeros([1, 5]))
    
    # Initial Design
    xt = param_choose(xt, get_dict_value(pypkpdInput, "design", "xt"), Exception("'xt' needs to be defined"))
    m = param_choose(m, get_dict_value(pypkpdInput, "design", "m"), None)
    x = param_choose(x, get_dict_value(pypkpdInput, "design", "x"), None)
    nx = param_choose(nx, get_dict_value(pypkpdInput, "design", "nx"), None)
    a = param_choose(a, get_dict_value(pypkpdInput, "design", "a"), None)
    groupsize = param_choose(groupsize, get_dict_value(pypkpdInput, "design", "groupsize"), Exception("'groupsize' needs to be defined"))
    ni = param_choose(ni, get_dict_value(pypkpdInput, "design", "ni"), None)
    model_switch = param_choose(model_switch, get_dict_value(pypkpdInput, "design", "model_switch"), None)
    
    # Design space
    maxni = param_choose(maxni, get_dict_value(pypkpdInput, "design_space", "maxni"), None)
    minni = param_choose(minni, get_dict_value(pypkpdInput, "design_space", "minni"), None)
    maxtotni = param_choose(maxtotni, get_dict_value(pypkpdInput, "design_space", "maxtotni"), None)
    mintotni = param_choose(mintotni, get_dict_value(pypkpdInput, "design_space", "mintotni"), None)
    maxgroupsize = param_choose(maxgroupsize, get_dict_value(pypkpdInput, "design_space", "maxgroupsize"), None)
    mingroupsize = param_choose(mingroupsize, get_dict_value(pypkpdInput, "design_space", "mingroupsize"), None)
    maxtotgroupsize = param_choose(maxtotgroupsize, get_dict_value(pypkpdInput, "design_space", "maxtotgroupsize"), None)
    mintotgroupsize = param_choose(mintotgroupsize, get_dict_value(pypkpdInput, "design_space", "mintotgroupsize"), None)
    maxxt = param_choose(maxxt, get_dict_value(pypkpdInput, "design_space", "maxxt"), None)
    minxt = param_choose(minxt, get_dict_value(pypkpdInput, "design_space", "minxt"), None)
    discrete_xt = param_choose(discrete_xt, get_dict_value(pypkpdInput, "design_space", "xt_space"), None)
    discrete_x = param_choose(discrete_x, get_dict_value(pypkpdInput, "design_space", "discrete_x"), None)
    maxa = param_choose(maxa, get_dict_value(pypkpdInput, "design_space", "maxa"), None)
    mina = param_choose(mina, get_dict_value(pypkpdInput, "design_space", "mina"), None)
    discrete_a = param_choose(discrete_a, get_dict_value(pypkpdInput, "design_space", "a_space"), None)
    bUseGrouped_xt = param_choose(bUseGrouped_xt, get_dict_value(pypkpdInput, "design_space", "bUseGrouped_xt"), False)
    G_xt = param_choose(G_xt, get_dict_value(pypkpdInput, "design_space", "G_xt"), None)
    bUseGrouped_a = param_choose(bUseGrouped_a, get_dict_value(pypkpdInput, "design_space", "bUseGrouped_a"), False)
    G_a = param_choose(G_a, get_dict_value(pypkpdInput, "design_space", "G_a"), None)
    bUseGrouped_x = param_choose(bUseGrouped_x, get_dict_value(pypkpdInput, "design_space", "bUseGrouped_x"), False)
    G_x = param_choose(G_x, get_dict_value(pypkpdInput, "design_space", "G_x"), None)

    # FIM calculation
    iFIMCalculationType = param_choose(iFIMCalculationType, get_dict_value(pypkpdInput, "settings", "iFIMCalculationType"), 1)
    iApproximationMethod = param_choose(iApproximationMethod, get_dict_value(pypkpdInput, "settings", "iApproximationMethod"), 0)
    iFOCENumInd = param_choose(iFOCENumInd, get_dict_value(pypkpdInput, "settings", "iFOCENumInd"), 1000)
    prior_fim = param_choose(prior_fim, get_dict_value(pypkpdInput, "settings", "prior_fim"), Matrix(np.zeros([0, 1])))
    strAutoCorrelationFile = param_choose(strAutoCorrelationFile, get_dict_value(pypkpdInput, "model", "auto_pointer"), "")

    # Criterion specification
    d_switch = param_choose(d_switch, get_dict_value(pypkpdInput, "settings", "prior_fim"), 1)
    ofv_calc_type = param_choose(ofv_calc_type, get_dict_value(pypkpdInput, "settings", "prior_fim"), 4)
    ds_index = param_choose(ds_index, get_dict_value(pypkpdInput, "parameters", "ds_index"), get_dict_value(pypkpdInput, "parameters", "ds_index"))
    strEDPenaltyFile = param_choose(strEDPenaltyFile, get_dict_value(pypkpdInput, "settings", "strEDPenaltyFile"), "")
    ofv_fun = param_choose(ofv_fun, get_dict_value(pypkpdInput, "settings", "ofv_fun"), None)
    
    # E-family Criterion options
    iEDCalculationType = param_choose(iEDCalculationType, get_dict_value(pypkpdInput, "settings", "iEDCalculationType"), 0)
    ED_samp_size = param_choose(ED_samp_size, get_dict_value(pypkpdInput, "settings", "ED_samp_size"), 45)
    bLHS = param_choose(bLHS, get_dict_value(pypkpdInput, "settings", "bLHS"), 1)
    strUserDistributionFile = param_choose(strUserDistributionFile, get_dict_value(pypkpdInput, "model", "user_distribution_pointer"), "")

    # Model parameters
    nbpop = param_choose(ds_index, get_dict_value(pypkpdInput, "parameters", "nbpop"), get_dict_value(pypkpdInput, "parameters", "nbpop"))
    NumRanEff = param_choose(NumRanEff, get_dict_value(pypkpdInput, "parameters", "NumRanEff"), get_dict_value(pypkpdInput, "parameters", "NumRanEff"))
    NumDocc = param_choose(NumDocc, get_dict_value(pypkpdInput, "parameters", "NumDocc"), get_dict_value(pypkpdInput, "parameters", "NumDocc"))
    NumOcc = param_choose(NumOcc, get_dict_value(pypkpdInput, "parameters", "NumOcc"), get_dict_value(pypkpdInput, "parameters", "NumOcc"))
    bpop = param_choose(bpop, get_dict_value(pypkpdInput, "parameters", "bpop"), Exception("bpop must be defined"))
    d = param_choose(d, get_dict_value(pypkpdInput, "parameters", "d"), None)
    covd = param_choose(covd, get_dict_value(pypkpdInput, "parameters", "covd"), get_dict_value(pypkpdInput, "parameters", "covd"))
    sigma = param_choose(sigma, get_dict_value(pypkpdInput, "parameters", "sigma"), get_dict_value(pypkpdInput, "parameters", "sigma"))
    docc = param_choose(docc, get_dict_value(pypkpdInput, "parameters", "docc"), Matrix(np.zeros([0, 3])))
    covdocc = param_choose(covdocc, get_dict_value(pypkpdInput, "parameters", "covdocc"), np.zeros([1, int(length(docc.get_data()[:, 1])*(length(docc.get_data()[:, 1])-1)/2)]))

    # Model parameters fixed or not
    notfixed_bpop = param_choose(notfixed_bpop, get_dict_value(pypkpdInput, "parameters", "notfixed_bpop"), get_dict_value(pypkpdInput, "parameters", "notfixed_bpop"))
    notfixed_d = param_choose(notfixed_d, get_dict_value(pypkpdInput, "parameters", "notfixed_d"), get_dict_value(pypkpdInput, "parameters", "notfixed_d"))
    notfixed_covd = param_choose(notfixed_covd, get_dict_value(pypkpdInput, "parameters", "notfixed_covd"), get_dict_value(pypkpdInput, "parameters", "notfixed_covd"))
    notfixed_docc = param_choose(notfixed_docc, get_dict_value(pypkpdInput, "parameters", "notfixed_docc"), get_dict_value(pypkpdInput, "parameters", "notfixed_docc"))
    notfixed_covdocc = param_choose(notfixed_covdocc, get_dict_value(pypkpdInput, "parameters", "notfixed_covdocc"), Matrix(np.zeros([1, length(covdocc)])))
    notfixed_sigma = param_choose(notfixed_sigma, get_dict_value(pypkpdInput, "parameters", "notfixed_sigma"), Matrix([1]*size(sigma)[1]))
    notfixed_covsigma = param_choose(notfixed_covsigma, get_dict_value(pypkpdInput, "parameters", "notfixed_covsigma"), Matrix(np.zeros([1, int(length(notfixed_sigma)*(length(notfixed_sigma)-1)/2)])))

    # Optimization algorithm choices
    bUseRandomSearch = param_choose(bUseRandomSearch, get_dict_value(pypkpdInput, "settings", "bUseRandomSearch"), True)
    bUseStochasticGradient = param_choose(bUseStochasticGradient, get_dict_value(pypkpdInput, "settings", "bUseStochasticGradient"), True)
    bUseLineSearch = param_choose(bUseLineSearch, get_dict_value(pypkpdInput, "settings", "bUseLineSearch"), True)
    bUseExchangeAlgorithm = param_choose(bUseExchangeAlgorithm, get_dict_value(pypkpdInput, "settings", "bUseExchangeAlgorithm"), False)
    bUseBFGSMinimizer = param_choose(bUseBFGSMinimizer, get_dict_value(pypkpdInput, "settings", "bUseBFGSMinimizer"), False)
    EACriteria = param_choose(EACriteria, get_dict_value(pypkpdInput, "settings", "EACriteria"), 1)
    strRunFile = param_choose(strRunFile, get_dict_value(pypkpdInput, "settings", "run_file_pointer"), "")

    # Labeling and file names
    pypkpd_version = param_choose(pypkpd_version, get_dict_value(pypkpdInput, "settings", "pypkpd_version"), None)
    modtit = param_choose(modtit, get_dict_value(pypkpdInput, "settings", "modtit"), "pypkpd model")
    output_file = param_choose(output_file, get_dict_value(pypkpdInput, "settings", "output_file"), "pypkpd_output_summary")
    output_function_file = param_choose(output_function_file, get_dict_value(pypkpdInput, "settings", "output_function_file"), "pypkpd_output")
    strIterationFileName = param_choose(strIterationFileName, get_dict_value(pypkpdInput, "settings", "strIterationFileName"), "pypkpd_current.py")
    
    # Misc options
    user_data = param_choose(user_data, get_dict_value(pypkpdInput, "settings", "user_data"), cell([0, 0]))
    ourzero = param_choose(ourzero, get_dict_value(pypkpdInput, "settings", "ourzero"), 1e-5)
    dSeed = param_choose(dSeed, get_dict_value(pypkpdInput, "settings", "dSeed"), None)
    line_opta = param_choose(line_opta, get_dict_value(pypkpdInput, "settings", "line_opta"), None)
    line_optx = param_choose(line_optx, get_dict_value(pypkpdInput, "settings", "line_optx"), None)
    bShowGraphs = param_choose(bShowGraphs, get_dict_value(pypkpdInput, "settings", "bShowGraphs"), False)
    use_logfile = param_choose(use_logfile, get_dict_value(pypkpdInput, "settings", "use_logfile"), False)
    m1_switch = param_choose(m1_switch, get_dict_value(pypkpdInput, "settings", "m1_switch"), 1)
    m2_switch = param_choose(m2_switch, get_dict_value(pypkpdInput, "settings", "m2_switch"), 1)
    hle_switch = param_choose(hle_switch, get_dict_value(pypkpdInput, "settings", "hle_switch"), 1)
    gradff_switch = param_choose(gradff_switch, get_dict_value(pypkpdInput, "settings", "gradff_switch"), 1)
    gradfg_switch = param_choose(gradfg_switch, get_dict_value(pypkpdInput, "settings", "gradfg_switch"), 1)
    grad_all_switch = param_choose(grad_all_switch, get_dict_value(pypkpdInput, "settings", "grad_all_switch"), 1)
    rsit_output = param_choose(rsit_output, get_dict_value(pypkpdInput, "settings", "rsit_output"), 5)
    sgit_output = param_choose(sgit_output, get_dict_value(pypkpdInput, "settings", "sgit_output"), 1)
    hm1 = param_choose(hm1, get_dict_value(pypkpdInput, "settings", "hm1"), 0.00001)
    hlf = param_choose(hlf, get_dict_value(pypkpdInput, "settings", "hlf"), 0.00001)
    hlg = param_choose(hlg, get_dict_value(pypkpdInput, "settings", "hlg"), 0.00001)
    hm2 = param_choose(hm2, get_dict_value(pypkpdInput, "settings", "hm2"), 0.00001)
    hgd = param_choose(hgd, get_dict_value(pypkpdInput, "settings", "hgd"), 0.00001)
    hle = param_choose(hle, get_dict_value(pypkpdInput, "settings", "hle"), 0.00001)
    AbsTol = param_choose(AbsTol, get_dict_value(pypkpdInput, "settings", "AbsTol"), 0.000001)
    RelTol = param_choose(RelTol, get_dict_value(pypkpdInput, "settings", "RelTol"), 0.000001)
    iDiffSolverMethod = param_choose(iDiffSolverMethod, get_dict_value(pypkpdInput, "settings", "iDiffSolverMethod"), None)
    bUseMemorySolver = param_choose(bUseMemorySolver, get_dict_value(pypkpdInput, "settings", "bUseMemorySolver"), False)
    rsit = param_choose(rsit, get_dict_value(pypkpdInput, "settings", "rsit"), 300)
    sgit = param_choose(sgit, get_dict_value(pypkpdInput, "settings", "sgit"), 150)
    intrsit = param_choose(intrsit, get_dict_value(pypkpdInput, "settings", "intrsit"), 250)
    intsgit = param_choose(intsgit, get_dict_value(pypkpdInput, "settings", "intsgit"), 50)
    maxrsnullit = param_choose(maxrsnullit, get_dict_value(pypkpdInput, "settings", "maxrsnullit"), 50)
    convergence_eps = param_choose(convergence_eps, get_dict_value(pypkpdInput, "settings", "convergence_eps"), 1e-08)
    rslxt = param_choose(rslxt, get_dict_value(pypkpdInput, "settings", "rslxt"), 10)
    rsla = param_choose(rsla, get_dict_value(pypkpdInput, "settings", "rsla"), 10)
    cfaxt = param_choose(cfaxt, get_dict_value(pypkpdInput, "settings", "cfaxt"), 0.001)
    cfaa = param_choose(cfaa, get_dict_value(pypkpdInput, "settings", "cfaa"), 0.001)
    bGreedyGroupOpt = param_choose(bGreedyGroupOpt, get_dict_value(pypkpdInput, "settings", "bGreedyGroupOpt"), False)
    EAStepSize = param_choose(EAStepSize, get_dict_value(pypkpdInput, "settings", "EAStepSize"), 0.01)
    EANumPoints = param_choose(EANumPoints, get_dict_value(pypkpdInput, "settings", "EANumPoints"), False)
    EAConvergenceCriteria = param_choose(EAConvergenceCriteria, get_dict_value(pypkpdInput, "settings", "EAConvergenceCriteria"), 1e-20)
    bEANoReplicates = param_choose(bEANoReplicates, get_dict_value(pypkpdInput, "settings", "bEANoReplicates"), False)
    BFGSProjectedGradientTol = param_choose(BFGSProjectedGradientTol, get_dict_value(pypkpdInput, "settings", "BFGSProjectedGradientTol"), 0.0001)
    BFGSTolerancef = param_choose(BFGSTolerancef, get_dict_value(pypkpdInput, "settings", "BFGSTolerancef"), 0.001)
    BFGSToleranceg = param_choose(BFGSToleranceg, get_dict_value(pypkpdInput, "settings", "BFGSToleranceg"), 0.9)
    BFGSTolerancex = param_choose(BFGSTolerancex, get_dict_value(pypkpdInput, "settings", "BFGSTolerancex"), 0.1)
    ED_diff_it = param_choose(ED_diff_it, get_dict_value(pypkpdInput, "settings", "ED_diff_it"), 30)
    ED_diff_percent = param_choose(ED_diff_percent, get_dict_value(pypkpdInput, "settings", "ED_diff_percent"), 10)
    line_search_it = param_choose(line_search_it, get_dict_value(pypkpdInput, "settings", "line_search_it"), 50)
    Doptim_iter = param_choose(Doptim_iter, get_dict_value(pypkpdInput, "settings", "Doptim_iter"), 1)

    # Parallel options for pypkpd
    iCompileOption = param_choose(iCompileOption, get_dict_value(pypkpdInput, "settings", "parallel", "iCompileOption"), -1)
    iUseParallelMethod = param_choose(iUseParallelMethod, get_dict_value(pypkpdInput, "settings", "parallel", "iUseParallelMethod"), 1)
    strExecuteName = param_choose(strExecuteName, get_dict_value(pypkpdInput, "settings", "parallel", "strExecuteName"), "calc_fim.exe")
    iNumProcesses = param_choose(iNumProcesses, get_dict_value(pypkpdInput, "settings", "parallel", "iNumProcesses"), 2)
    iNumChunkDesignEvals = param_choose(iNumChunkDesignEvals, get_dict_value(pypkpdInput, "settings", "parallel", "iNumChunkDesignEvals"), -2)
    Mat_Out_Pre = param_choose(Mat_Out_Pre, get_dict_value(pypkpdInput, "settings", "parallel", "strMatFileOutputPrefix"), "parallel_output")
    strExtraRunOptions = param_choose(strExtraRunOptions, get_dict_value(pypkpdInput, "settings", "parallel", "strExtraRunOptions"), "")
    dPollResultTime = param_choose(dPollResultTime, get_dict_value(pypkpdInput, "settings", "parallel", "dPollResultTime"), 0.1)
    strFunctionInputName = param_choose(strFunctionInputName, get_dict_value(pypkpdInput, "settings", "parallel", "strFunctionInputName"), "function_input")
    bParallelRS = param_choose(bParallelRS, get_dict_value(pypkpdInput, "settings", "parallel", "bParallelRS"), False)
    bParallelSG = param_choose(bParallelSG, get_dict_value(pypkpdInput, "settings", "parallel", "bParallelSG"), False)
    bParallelMFEA = param_choose(bParallelMFEA, get_dict_value(pypkpdInput, "settings", "parallel", "bParallelMFEA"), False)
    bParallelLS = param_choose(bParallelLS, get_dict_value(pypkpdInput, "settings", "parallel", "bParallelLS"), False)

#####-------------- main part --------------#####
    pypkpd_db = {}

    pypkpd_db["settings"] = {}
    pypkpd_db["settings"]["pypkpd_version"] = pypkpd_version
    
    if BFGSConvergenceCriteriaMinStep is None:
        BFGSConvergenceCriteriaMinStep = param_choose(BFGSConvergenceCriteriaMinStep, get_dict_value(pypkpdInput, "settings", "BFGSConvergenceCriteriaMinStep"), 1e-08)
    if MCC_Dep is None:
        MCC_Dep = param_choose(MCC_Dep, get_dict_value(pypkpdInput, "settings", "parallel", "strAdditionalMCCCompilerDependencies"), "")

    pypkpd_db["model"] = {}
    pypkpd_db["model"]["user_distribution_pointer"] = ""

    # if strUserDistributionFile is not None:
    #     if str(strUserDistributionFile) != "":
    #         pypkpd_db["model"]["user_distribution_pointer"] = strUserDistributionFile
    # else:
    #     exec(open(strUserDistributionFile).read())
    #     returnArgs = fileparts(strUserDistributionFile)
    #     strUserDistFilePath = list(returnArgs.values())[0]
    #     strUserDistFilename = list(returnArgs.values())[1]
    #     pypkpd_db["model"]["user_distribution_pointer"] = strUserDistFilename

    
    if np.any(np.array(size(x)) == 0): # should be removed
        x = None
        G_x = None
        discrete_x = None

    if np.any(np.array(size(a)) == 0): # should be removed
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
    if get_dict_value(design_space, "use_grouped_a") is not None:
        design_space["bUseGrouped_a"] = get_dict_value(design_space, "use_grouped_a")
        design_space["use_grouped_a"] = None

    if get_dict_value(design_space, "use_grouped_x") is not None:
        design_space["bUseGrouped_x"] = get_dict_value(design_space, "use_grouped_x")
        design_space["use_grouped_x"] = None

    if get_dict_value(design_space, "use_grouped_xt") is not None:
        design_space["bUseGrouped_xt"] = get_dict_value(design_space, "use_grouped_xt")
        design_space["use_grouped_xt"] = None

    if get_dict_value(design_space, "grouped_a") is not None:
        design_space["G_a"] = get_dict_value(design_space, "grouped_a")
        design_space["grouped_a"] = None

    if get_dict_value(design_space, "grouped_x") is not None:
        design_space["G_x"] = get_dict_value(design_space, "grouped_x")
        design_space["grouped_x"] = None

    if get_dict_value(design_space, "grouped_xt") is not None:
        design_space["G_xt"] = get_dict_value(design_space, "grouped_xt")
        design_space["grouped_xt"] = None

    if get_dict_value(design_space, "x_space") is not None:
        design_space["discrete_x"] = get_dict_value(design_space, "x_space")
        design_space["x_space"] = None

    # should be removed
    if get_dict_value(design, "x") is None:
        design["x"] = Matrix(np.zeros([get_dict_value(design, "m").get_value(), 0]))
        design_space["G_x"] = design["x"]
        design_space["bUseGrouped_x"] = False
        design_space["discrete_x"] = cell([get_dict_value(design, "m").get_value(), 0], fill_value=1)

    # should be removed
    if get_dict_value(design, "a") is None:
        design["a"] = Matrix(np.zeros([get_dict_value(design, "m").get_value(), 0]))
        design_space["G_a"] = design["a"]
        design_space["bUseGrouped_a"] = False
        design_space["mina"] = design["a"]
        design_space["maxa"] = design["a"]

    pypkpd_db["design"] = design
    pypkpd_db["design_space"] = design_space

    pypkpd_db["settings"]["bLHS"] = bLHS

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
    pypkpd_db["settings"]["solved_solutions"] = cell([0, 0])
    pypkpd_db["settings"]["maxtime"] = max(get_dict_value(pypkpd_db, "design_space", "maxxt").get_data()) + get_dict_value(pypkpd_db, "settings", "hgd")

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
        pypkpd_db["model"]["fg_pointer"] = fg_fun
    elif fg_file is not None:
        pypkpd_db["model"]["fg_pointer"] = fg_file
    # omitted
    # else:
    #     try:
    #         fg_file
    #     except NameError:
    #         # exec(open(fg_file).read())
    #         returnArgs = fileparts(fg_file)
    #         strfgModelFilePath = list(returnArgs.values())[0]
    #         strfgModelFilename = list(returnArgs.values())[1]
    #         # if (~strcmp(strfgModelFilePath,''))
    #         # cd(strfgModelFilePath);
    #         # end
    #         pypkpd_db["model"]["fg_pointer"] = strfgModelFilename
    #     else:
    #         pypkpd_db["model"]["fg_pointer"] = fg_file

    pypkpd_db["settings"]["ed_penalty_pointer"] = np.zeros([1, 0])
    
    # omitted
    # if strEDPenaltyFile is not None:
    #     if str(strEDPenaltyFile) != "":
    #         pypkpd_db["settings"]["ed_penalty_pointer"] = strEDPenaltyFile
    
    if (ofv_fun is None) or callable(ofv_fun):
        pypkpd_db["settings"]["ofv_fun"] = ofv_fun
    # omitted

    pypkpd_db["model"]["auto_pointer"] = np.zeros([1, 0])
    
    # omitted
    # if strAutoCorrelationFile != "":
    #     if str(strAutoCorrelationFile) != "":
    #         if strAutoCorrelationFile is not None:
    #             pypkpd_db["model"]["auto_pointer"] = strAutoCorrelationFile
    #         else:
    #             exec(open(pypkpdInput["strAutoCorrelationFile"]).read())
    #             returnArgs = fileparts(pypkpdInput["strAutoCorrelationFile"])
    #             strAutoCorrelationFilePath = list(returnArgs.values())[0]
    #             strAutoCorrelationFilename = list(returnArgs.values())[1]
    #         # if (~strcmp(strAutoCorrelationFilePath,''))
    #         # cd(strAutoCorrelationFilePath);
    #         # end
    #             pypkpd_db["model"]["auto_pointer"] = strAutoCorrelationFilename

    if callable(ff_fun):
        pypkpd_db["model"]["ff_pointer"] = ff_fun
    elif type(ff_fun) is str:
        pypkpd_db["model"]["ff_pointer"] = ff_fun
    elif ff_file is not None:
        pypkpd_db["model"]["ff_pointer"] = ff_file
    # omitted
    # else:
    #     try:
    #         ff_file
    #     except NameError:
    #         # exec(open(ff_file).read())
    #         returnArgs = fileparts(ff_file)
    #         strffModelFilePath = list(returnArgs.values())[0]
    #         strffModelFilename = list(returnArgs.values())[1]
    #         # if (~strcmp(strffModelFilePath,''))
    #         # cd(strffModelFilePath);
    #         # end
    #         pypkpd_db["model"]["ff_pointer"] = strffModelFilename
    #     else:
    #         pypkpd_db["model"]["ff_pointer"] = ff_file

    # omitted
    # Check if there is any sub models defined
    # if "SubModels" in pypkpdInput.keys():
    #     i = 1
    #     for elem_key in pypkpdInput["SubModels"]:
    #         if elem_key == ("ff_file%d" % i):
    #             # ok<np.nanSGU>
    #             exec(eval('pypkpdInput["SubModels"]["ff_file"]%d' % i))
    #             returnArgs = fileparts(eval('pypkpdInput["SubModels"]["ff_file"]%d' % i))  # ok<np.nanSGU>
    #             strffModelFilePath = list(returnArgs.values())[0]
    #             strffModelFilename = list(returnArgs.values())[1]
    #             # if (~strcmp(strffModelFilePath,''))
    #             # cd(strffModelFilePath);
    #             # end
    #             pypkpd_db["model"]["subffPointers"]["".join(["ff_pointer", str(i)])] = strffModelFilename
    #             i = i + 1

    if callable(fError_fun):
        pypkpd_db["model"]["ferror_pointer"] = fError_fun
    elif type(fError_fun) is str:
        pypkpd_db["model"]["ferror_pointer"] = fError_fun
    elif fError_file is not None:
        pypkpd_db["model"]["fError_pointer"] = fError_file
    # omitted
    # else:
    #     try:
    #         fError_file
    #     except NameError:
    #         # exec(open(fError_file).read())
    #         returnArgs = fileparts(fError_file)
    #         strErrorModelFilePath = list(returnArgs.values())[0]
    #         strErrorModelFilename = list(returnArgs.values())[1]
    #         # if (~strcmp(strErrorModelFilePath,''))
    #         # cd(strErrorModelFilePath);
    #         # end
    #         pypkpd_db["model"]["ferror_pointer"] = strErrorModelFilename
    #     else:
    #         pypkpd_db["model"]["ferror_pointer"] = fError_file


# ==================================
# Initialize the randomization
# ==================================
    # omitted
    # if dSeed is not None:
    #     if dSeed == -1:
    #         pypkpd_db["settings"]["dSeed"] = datetime.datetime.now()
    #     else:
    #         pypkpd_db["settings"]["dSeed"] = dSeed
    #     random.seed(pypkpd_db["settings"]["dSeed"])

    pypkpd_db["parameters"]["nbpop"] = pypkpd_choose(nbpop, find_largest_index(get_dict_value(pypkpd_db, "model", "fg_pointer"), "bpop"))
    pypkpd_db["parameters"]["NumRanEff"] = pypkpd_choose(NumRanEff, find_largest_index(get_dict_value(pypkpd_db, "model", "fg_pointer"), "b"))
    pypkpd_db["parameters"]["NumDocc"] = pypkpd_choose(NumDocc, find_largest_index(get_dict_value(pypkpd_db, "model", "fg_pointer"), "bocc", mat=True, mat_row=True))
    pypkpd_db["parameters"]["NumOcc"] = pypkpd_choose(NumOcc, find_largest_index(get_dict_value(pypkpd_db, "model", "fg_pointer"), "bocc", mat=True, mat_row=False))
    
    tmp_ng = call_func(get_dict_value(pypkpd_db, "model", "fg_pointer"), 0, 0, 0, 0, np.zeros([get_dict_value(pypkpd_db, 'parameters', 'NumDocc'), get_dict_value(pypkpd_db, 'parameters', 'NumOcc')]))
    if type(tmp_ng) is Matrix:
        pypkpd_db["parameters"]["ng"] = tmp_ng.get_size()
    else:
        raise Exception("ng needs to be defined.")
    
    pypkpd_db["parameters"]["notfixed_docc"] = pypkpd_choose(notfixed_docc, Matrix(np.ones([1, get_dict_value(pypkpd_db, "parameters", "NumDocc")])))
    pypkpd_db["parameters"]["notfixed_d"] = pypkpd_choose(notfixed_d, Matrix(np.ones([1, get_dict_value(pypkpd_db, "parameters", "NumRanEff")])))
    pypkpd_db["parameters"]["notfixed_bpop"] = pypkpd_choose(notfixed_bpop, Matrix(np.ones([1, get_dict_value(pypkpd_db, "parameters", "nbpop")])))


    # reorder named values
    fg_names = call_func(get_dict_value(pypkpd_db, "model", "fg_pointer"), 0, 0, 0, 0, np.ones([get_dict_value(pypkpd_db, 'parameters', 'NumDocc'), get_dict_value(pypkpd_db, 'parameters', 'NumOcc')])).get_datanam()
    
    pypkpd_db["parameters"]["notfixed_bpop"] = reorder_vec(get_dict_value(pypkpd_db, "parameters", "notfixed_bpop"), fg_names)

    pypkpd_db["parameters"]["notfixed_d"] = reorder_vec(get_dict_value(pypkpd_db, "parameters", "notfixed_d"), fg_names)

    # we have just the parameter values not the uncertainty
    if d.get_shape(0) == 1 and d.get_shape(1) == get_dict_value(pypkpd_db, "parameters", "NumRanEff"):
        d_descr = np.zeros([get_dict_value(pypkpd_db, "parameters", "NumRanEff"), 3])

        # reorder named values
        d = reorder_vec(d, fg_names)
        
        d_descr[:, 1] = d.get_data()
        d_descr[:, 0] = np.zeros(d_descr.shape[0])
        d_descr[:, 2] = np.zeros(d_descr.shape[0])
        d_descr = Matrix(d_descr, axisnam=[d.get_datanam(), None])
        d = d_descr

    # we have just the parameter values not the uncertainty
    if bpop.get_shape(0) == 1 and bpop.get_shape(1) == get_dict_value(pypkpd_db, "parameters", "nbpop"):
        bpop_descr = np.zeros([get_dict_value(pypkpd_db, "parameters", "nbpop"), 3])

        # reorder named values
        bpop = reorder_vec(bpop, fg_names)
    
        bpop_descr[:, 1] = bpop.get_data()
        bpop_descr[:, 0] = np.zeros(bpop_descr.shape[0])
        bpop_descr[:, 2] = np.zeros(bpop_descr.shape[0])
        bpop_descr = Matrix(bpop_descr, axisnam=[bpop.get_datanam(), None])
        bpop = bpop_descr

    # we have just the diagonal parameter values
    if type(sigma) is not Matrix and size(sigma)[0] == 1:
        sigma_tmp = Matrix(np.diag([sigma] * size(sigma)[1]))
        sigma_tmp.set_axisnam([names(sigma), None])
        sigma = sigma_tmp

    covd = pypkpd_choose(covd, np.zeros([1, int(get_dict_value(pypkpd_db, "parameters", "NumRanEff")*(get_dict_value(pypkpd_db, "parameters", "NumRanEff")-1)/2)]))
    pypkpd_db["parameters"]["covd"] = covd

    tmp = np.ones(1, length(covd))
    tmp[covd.get_data() == 0] = np.zeros(np.sum(covd.get_data() == 0))
    pypkpd_db["parameters"]["notfixed_covd"] = pypkpd_choose(notfixed_covd, tmp)

    # ==================================
    # Sample the individual eta's for FOCE and FOCEI
    # ==================================

    if get_dict_value(pypkpd_db, "settings", "iApproximationMethod") != 0 and get_dict_value(pypkpd_db, "settings", "iApproximationMethod") != 3:

        iMaxCorrIndNeeded = 100
        # bzeros = np.zeros([get_dict_value(pypkpd_db, "parameters", "NumRanEff"), 1])
        # bones = np.ones([get_dict_value(pypkpd_db, "parameters", "NumRanEff"), 1])
        # bocczeros = np.zeros([get_dict_value(pypkpd_db, "parameters", "NumDocc"), 1])
        # boccones = np.ones([get_dict_value(pypkpd_db, "parameters", "NumDocc"), 1])

        pypkpd_db["parameters"]["b_global"] = np.zeros([get_dict_value(pypkpd_db, "parameters", "NumRanEff"), max(get_dict_value(pypkpd_db, "settings", "iFOCENumInd"), iMaxCorrIndNeeded)])

        fulld = getfulld(d.get_partial_matrix([[None, None], [1, 1]]), get_dict_value(pypkpd_db, "parameters", "covd"))
        fulldocc = getfulld(Matrix(docc.get_partial_matrix([[None, None], [1, 1]])), get_dict_value(pypkpd_db, "parameters", "covdocc"))

        pypkpd_db["parameters"]["bocc_global"] = cell([get_dict_value(pypkpd_db, "settings", "iFOCENumInd"), 1], fill_value=0)

        if get_dict_value(pypkpd_db, "settings", "d_switch") is True:
            pypkpd_db["parameters"]["b_global"] = Matrix(np.transpose(np.random.multivariate_normal(max(get_dict_value(pypkpd_db, "settings", "iFOCENumInd"), iMaxCorrIndNeeded), sigma=fulld)))
            for i in range(0, get_dict_value(pypkpd_db, "settings", "iFOCENumInd")):
                pypkpd_db["parameters"]["bocc_global"][i] = np.zeros([docc.get_shape(0), get_dict_value(pypkpd_db, "parameters", "NumOcc")])
                if get_dict_value(pypkpd_db, "parameters", "NumOcc") != 0:
                    pypkpd_db["parameters"]["bocc_global"][i] = np.transpose(np.random.multivariate_normal(get_dict_value(pypkpd_db, "parameters", "NumOcc"), sigma=fulldocc))

        else:
            d_dist = pargen(d, get_dict_value(pypkpd_db, "model", "user_distribution_pointer"), max(get_dict_value(pypkpd_db, "settings", "iFOCENumInd"), iMaxCorrIndNeeded), get_dict_value(pypkpd_db, "settings", "bLHS"), np.zeros([1, 0]), pypkpd_db)
            docc_dist = pargen(docc, get_dict_value(pypkpd_db, "model", "user_distribution_pointer"), get_dict_value(pypkpd_db, "settings", "iFOCENumInd"), get_dict_value(pypkpd_db, "settings", "bLHS"), np.zeros([1, 0]), pypkpd_db)

            ## set one data with type of Matrix
            if d_dist.get_size() != 0:
                for i in range(0, max(get_dict_value(pypkpd_db, "settings", "iFOCENumInd"), iMaxCorrIndNeeded)):
                    for j in range(0, get_dict_value(pypkpd_db, "parameters", "b_global").get_shape(0)):
                        get_dict_value(pypkpd_db, "parameters", "b_global").set_one_data(Matrix(np.transpose(np.random.multivariate_normal(1, sigma=getfulld(d_dist[i, :], get_dict_value(pypkpd_db, "parameters", "covd"))))), index=[j, i])                        

            if docc_dist.get_size() != 0:
                for i in range(0, get_dict_value(pypkpd_db, "settings", "iFOCENumInd")):
                    pypkpd_db["parameters"]["bocc_global"].set_one__data(Matrix(np.transpose(np.random.multivariate_normal(get_dict_value(pypkpd_db, "parameters", "NumOcc"), sigma=getfulld(docc_dist[i, :], get_dict_value(pypkpd_db, "parameters", "covdocc"))))), index=[0, i])
    else:
        pypkpd_db["parameters"]["bocc_global"] = np.zeros([get_dict_value(pypkpd_db, "parameters", "NumRanEff"), 1])
        pypkpd_db["parameters"]["bocc_global"] = cell([1, 1], fill_value=1)
        pypkpd_db["parameters"]["bocc_global"].set_one_data(np.zeros(size(docc)[0], get_dict_value(pypkpd_db, "parameters", "NumOcc")), index=[0,1])
        
        pypkpd_db["settings"]["iFOCENumInd"] = 1

    pypkpd_db["settings"]["modtit"] = modtit
    pypkpd_db["settings"]["exptit"] = ('%s_exp["mat"]', modtit)  # experiment settings title
    pypkpd_db["settings"]["opttit"] = ('%s_opt["mat"]', modtit)  # optimization settings title
    pypkpd_db["settings"]["bShowGraphs"] = bShowGraphs

    pypkpd_db["settings"]["use_logfile"] = use_logfile
    pypkpd_db["settings"]["output_file"] = output_file
    pypkpd_db["settings"]["output_function_file"] = output_function_file

    pypkpd_db["settings"]["optsw"] = optsw

    line_opta = pypkpd_choose(line_opta, np.ones([1, get_dict_value(pypkpd_db, "design", "a").get_shape(1)]))
    if test_mat_size(np.array([1, get_dict_value(pypkpd_db, "design", "a").get_shape(1)]), line_opta, "line_opta"):
        pypkpd_db["settings"]["line_opta"] = Matrix(line_opta)

    line_optx = pypkpd_choose(line_optx, np.ones([1, get_dict_value(pypkpd_db, "design", "x").get_shape(1)]))
    if test_mat_size(np.array([1, get_dict_value(pypkpd_db, "design", "x").get_shape(1)]), line_optx, "line_optx"):
        pypkpd_db["settings"]["line_optx"] = Matrix(line_optx)

    pypkpd_db["parameters"]["bpop"] = bpop
    pypkpd_db["parameters"]["d"] = d
    pypkpd_db["parameters"]["covd"] = covd
    pypkpd_db["parameters"]["sigma"] = sigma
    pypkpd_db["parameters"]["docc"] = docc
    pypkpd_db["parameters"]["covdocc"] = covdocc

    pypkpd_db["settings"]["m1_switch"] = m1_switch
    pypkpd_db["settings"]["m2_switch"] = m2_switch
    pypkpd_db["settings"]["hle_switch"] = hle_switch
    pypkpd_db["settings"]["gradff_switch"] = gradff_switch
    pypkpd_db["settings"]["gradfg_switch"] = gradfg_switch
    pypkpd_db["settings"]["grad_all_switch"] = grad_all_switch

    pypkpd_db["settings"]["prior_fim"] = prior_fim

    pypkpd_db["parameters"]["notfixed_sigma"] = notfixed_sigma

    pypkpd_db["parameters"]["ds_index"] = ds_index

    # create ds_index if not already done
    if get_dict_value(pypkpd_db, "parameters", "ds_index") is None:
        unfixed_params = get_unfixed_params(pypkpd_db)
        pypkpd_db["parameters"]["ds_index"] = Matrix(np.transpose(np.repeat([0], unfixed_params.get_one_data(name="all").size)))
        pypkpd_db["parameters"]["ds_index"].data[(unfixed_params.get_one_data(name="bpop").size+1):pypkpd_db["parameters"]["ds_index"].get_size()] = 1
    else:
        if type(pypkpd_db["parameters"]["ds_index"]) is not Matrix:
            pypkpd_db["parameters"]["ds_index"] = Matrix(pypkpd_db["parameters"]["ds_index"], shape=[1, pypkpd_db["parameters"]["ds_index"].size])

    pypkpd_db["settings"]["strIterationFileName"] = strIterationFileName
    pypkpd_db["settings"]["user_data"] = user_data
    pypkpd_db["settings"]["bUseSecondOrder"] = False
    pypkpd_db["settings"]["bCalculateEBE"] = False
    pypkpd_db["settings"]["bGreedyGroupOpt"] = bGreedyGroupOpt
    pypkpd_db["settings"]["bEANoReplicates"] = bEANoReplicates

    if len(strRunFile) != 0:
        if strRunFile == " ":
            pypkpd_db["settings"]["run_file_pointer"] = np.zeros(1, 0)
    #     else:
    #         if strRunFile is not None:
    #             pypkpd_db["settings"]["run_file_pointer"] = strRunFile
    #         else:
    #             exec(open(pypkpdInput["strRunFile"]).read())
    #             returnArgs = fileparts(pypkpdInput["strRunFile"])
    #             strRunFilePath = list(returnArgs.values())[0]
    #             strRunFilename = list(returnArgs.values())[1]
    #             pypkpd_db["settings"]["run_file_pointer"] = strRunFilename

    # pypkpd_db["settings"]["Engine"] = {"Type": 1, "Version": version["version.string"]}
    
    pypkpd_db = convert_variables(pypkpd_db)  # need to transform here

    param_val = get_all_params(pypkpd_db)
    tmp_names = param_val.get_datanam()
    eval('%s.val = param.val["%s"]' % (tmp_names, tmp_names))
    
    d_val = d_val # for package check
    covd_val = covd_val
    docc_val = docc_val
    covdocc_val = covdocc_val
    bpop_val = bpop_val
    d_full = getfulld(d_val, covd_val)
    docc_full = getfulld(docc_val, covdocc_val)
    sigma_full = pypkpd_db["parameters"]["sigma"]
    
    pypkpd_db["parameters"]["param_pt_val"]["bpop"] = bpop_val
    pypkpd_db["parameters"]["param_pt_val"]["d"] = d_full
    pypkpd_db["parameters"]["param_pt_val"]["docc"] = docc_full
    pypkpd_db["parameters"]["param_pt_val"]["sigma"] = sigma_full

    retargs = fileparts(pypkpd_db["settings"]["output_file"])
    pypkpd_db["settings"]["strOutputFilePath"] = list(retargs.values())[0]
    pypkpd_db["settings"]["strOutputFileName"] = list(retargs.values())[1]
    pypkpd_db["settings"]["strOutputFileExtension"] = list(retargs.values())[2]

    return pypkpd_db

def somestring(**kwargs):
    return ", ".join(f"{key}={value}" for key, value in kwargs.items())