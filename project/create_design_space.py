"""
#' Create design variables and a design space for a full description of an optimization problem.
#' 
#' \code{create_design_space} takes an initial design and arguments for a design space and 
#' creates a design and design space for design optimization.
#' Checks the sizes of supplied design space variables and 
#' changes them to sizes that  make sense if there are inconsistencies.  
#' Function arguments can use shorthand notation (single values, vectors, lists of vectors and 
#' list of list) or matricies.
#' Returns a list of matricies compatible with PopED.
#' 
#' If a value (or a vector or a list of values) is supplied that corresponds to only one group and the design has
#' multiple groups then all groups will have the same value(s). If a matrix is expected then a list of lists can be supplied 
#' instead, each list corresponding to a group.   
#' 
#' @param design  The output from a call to \code{\link{create_design}}.
#' @param maxni Vector defining the maximum number of samples per group. 
#' @param minni Vector defining the minimum number of samples per group. 
#' @param maxtotni Number defining the maximum number of samples allowed in the experiment. 
#' @param mintotni Number defining the minimum number of samples allowed in the experiment.  
#' @param maxgroupsize Vector defining the maximum size of the different groups (maximum number of individuals in each group)
#' @param mingroupsize Vector defining the minimum size of the different groups (minimum num individuals in each group) 
#' @param maxtotgroupsize The total maximal groupsize over all groups
#' @param mintotgroupsize The total minimal groupsize over all groups
#' @param maxxt Matrix or single value defining the maximum value for each xt sample.  If a single value is 
#' supplied then all xt values are given the same maximum value.
#' @param minxt Matrix or single value defining the minimum value for each xt sample.  If a single value is 
#' supplied then all xt values are given the same minimum value
#' @param x_space Cell array \code{\link{cell}} defining the discrete variables for each x value.
#' @param xt_space Cell array \code{\link{cell}} defining the discrete variables allowed for each xt value.
#'   Can also be a vector of values \code{c(1:10)} (same values allowed for all xt), or a list of lists 
#'  \code{list(1:10, 2:23, 4:6)} (one for each value in xt in row major order or just for one row in xt, 
#'  and all other rows will be duplicated).
#' @param a_space Cell array \code{\link{cell}} defining the discrete variables allowed for each a value.
#'   Can also be a list of values \code{list(1:10)} (same values allowed for all a), or a list of lists 
#'  \code{list(1:10, 2:23, 4:6)} (one for each value in a).
#' @param maxa Vector defining the maximum value for each covariate. IF a single value is supplied then
#'  all a values are given the same maximum value
#' @param mina Vector defining the minimum value for each covariate. IF a single value is supplied then
#'  all a values are given the same minimum value
#' @param use_grouped_xt Group sampling times between groups so that each group has the same values (\code{TRUE} or \code{FALSE}).
#' @param grouped_xt Matrix defining the grouping of sample points. Matching integers mean that the points are matched.  
#' Allows for finer control than \code{use_grouped_xt}
#' @param use_grouped_a Group continuous design variables between groups so that each group has the same values (\code{TRUE} or \code{FALSE}).
#' @param grouped_a Matrix defining the grouping of continuous design variables. Matching integers mean that the values are matched.  
#' Allows for finer control than \code{use_grouped_a}.
#' @param use_grouped_x Group discrete design variables between groups so that each group has the same values (\code{TRUE} or \code{FALSE}).
#' @param grouped_x Matrix defining the grouping of discrete design variables. Matching integers mean that the values are matched.  
#' Allows for finer control than \code{use_grouped_x}.
#' @param our_zero Value to interpret as zero in design.
#' 
#' 
#' @family poped_input
#' 
#' @example tests/testthat/examples_fcn_doc/examples_create_design_space.R
#' 
#' @export
#' 
# @importFrom dplyr rbind_all

Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
from numpy.core.fromnumeric import shape
import pandas as pd
import inspect
from project.size import size
from project.test_mat_size import test_mat_size
from project.ones import ones
from project.cell import cell


def create_design_space(design,
						maxni=None,
						minni=None,
						maxtotni=None,
						mintotni=None,
						maxgroupsize=None,
						mingroupsize=None,
						maxtotgroupsize=None,
						mintotgroupsize=None,
						maxxt=None,
						minxt=None,
						xt_space=None,
						maxa=None,
						mina=None,
						a_space=None,
						x_space=None,
						use_group_xt=False,
						grouped_xt=None,
						use_grouped_a=False,
						grouped_a=None,
						use_grouped_x=False,
						grouped_x=None,
						our_zero=None):

	called_args = inspect.getargspec(create_design_space).args

	# assign defaults if not supplied
	if maxni is None:
		maxni = design["ni"]
	if minni is None:
		minni = design["ni"]
	if maxgroupsize is None:
		maxgroupsize = design["groupsize"]
	if mingroupsize is None:
		mingroupsize = design["groupsize"]
	maxxt_imputed = False
	if maxxt is None:
		maxxt = design["xt"]
		maxxt_imputed = True
	minxt_imputed = False
	if minxt is None:
		minxt = design["xt"]
		minxt_imputed = True
	maxa_imputed = False
	if maxa is None:
		maxa = design["a"]
		maxa_imputed = True
	mina_imputed = False
	if mina is None:
		mina = design["a"]
		mina_imputed = True
	
	design_space = {}
	design_new = design

	# rules:
	# 1. set default if not already defined
	# 2. read in value and translate to correct format
	# 3. check that size of object is correct
	# 4. add row and column names
	# 4. check that max is greater than min
	# 5. check that design value is wihin range of design_space

	# maxni
	if size(maxni, 1) == 1 and design["m"] != 1:
		maxni = np.array([maxni] * design["m"]).reshape([design["m"], 1])
	if type(maxni) is not np.ndarray:
		maxni = np.array(maxni)
	if test_mat_size(np.array([design["m"], 1]), maxni, "maxni") == 1:
		maxni = pd.DataFrame(maxni,
							 index=["grp_"+str(i) for i in range(1, design["m"]+1)],
							 columns="n_obs"*maxni.shape[1])

	# minni
	if size(minni, 1) == 1 and design["m"] != 1:
		minni = np.array([minni] * design["m"]).reshape([design["m"], 1])
	if type(minni) is not np.ndarray:
		minni = np.array(minni)
	if test_mat_size(np.array([design["m"], 1]), minni, "minni") == 1:
		minni = pd.DataFrame(minni,
							 index=["grp_"+str(i) for i in range(1, design["m"]+1)],
							 columns="n_obs"*minni.shape[1])
	
	# make sure min is smaller than max
	ret = comp_max_min(maxni, minni, called_args)
	maxni = ret["max_val"]
	minni = ret["min_val"]

	# check ni given max and min
	if any(design["ni"] < minni):
		raise Exception("ni is less than minni")
	if any(design["ni"] > maxni):
		raise Exception("ni is greater than maxni")
	
	# maxtotni and mintotni
	if maxtotni is None:
		maxtotni = np.sum(maxni)
	if mintotni is None:
		mintotni = np.sum(minni)
	test_mat_size(np.array([1, 1]), maxtotni, "maxtotni")
	test_mat_size(np.array([1, 1]), mintotni, "mintotni")
	ret = comp_max_min(maxtotni, mintotni, called_args)
	maxtotni = ret["max_val"]
	mintotni = ret["min_val"]
	if any(np.sum(design["ni"]) < mintotni):
		raise Exception("sum of ni is less than mintotni")
	if any(np.sum(design["ni"]) > maxtotni):
		raise Exception("sum of ni is greater than maxtotni")

	# update xt and model_switch given maxni
	if np.amax(maxni) > size(design["xt"], 2):

		# xt has to increase
		xt_full = ones(design["m"], np.amax(maxni)) * np.nan
		xt_full[0:design["m"], 0:size(design["xt"], 2)] = design["xt"]
		xt_full = pd.DataFrame(xt_full,
							   index=["grp_"+str(i) for i in range(1, design["m"]+1)],
							   columns=["obs_"+str(i) for i in range(1, size(xt_full, 2)+1)])
		design["xt"] = xt_full
		design_new["xt"] = design["xt"]

		# model switch has to increase
		model_switch_full = ones(design["m"], np.amax(maxni)) * np.nan
		model_switch_full[0:design["m"], 0:size(design["model_switch"], 2)] = design["model_switch"]
		model_switch_full = pd.DataFrame(model_switch_full,
							   index=["grp_"+str(i) for i in range(1, design["m"]+1)],
							   columns=["obs_"+str(i) for i in range(1, size(model_switch_full, 2)+1)])
		design["model_switch"] = model_switch_full
		for i in range(1, design["model_switch"].shape[0]+1):
			x_tmp = design["model_switch"][i, :]
			_, idx = np.unique(x_tmp[~np.isnan(x_tmp)], return_index=True) # remove duplicated and nan values but keep order
			if len([x_tmp[~np.isnan(x_tmp)][index] for index in sorted(idx)]) == 1:
				x_tmp[np.isnan(x_tmp)] = [x_tmp[~np.isnan(x_tmp)][index] for index in sorted(idx)]
			else:
				raise Exception("Unable to determine the model_switch values needed for group " + str(i)
								+ "\n Please supply them as input.")
			design["model_switch"][i, :] = x_tmp
		design_new["model_switch"] = design["model_switch"]

	# maxgroupsize
	if size(maxgroupsize, 1) == 1 and design["m"] != 1:
		...
		# maxgroupsize = np.array()
	
	# mingroupsize


	# make sure min is less than max
	ret = comp_max_min(maxgroupsize, mingroupsize, called_args)
	maxgroupsize = ret["max_val"]
	mingroupsize = ret["min_val"]

	# check given max and min
	if any(design["groupsize"] < mingroupsize):
		raise Exception("groupsize is less than mingroupsize")
	if any(design["groupsize"] > maxgroupsize):
		raise Exception("groupsize is greater than maxgroupsize")

	# maxtotgroupsize
	if maxtotgroupsize is None:
		maxtotgroupsize = np.sum(design["groupsize"])
	
	# mintotgroupsize
	if mintotgroupsize is None:
		mintotgroupsize = np.sum(mingroupsize)

	# make sure min is less than max
	ret = comp_max_min(maxtotgroupsize, mintotgroupsize, called_args)
	maxtotgroupsize = ret["max_val"]
	mintotgroupsize = ret["min_val"]

	# check given max and min
	if any(np.sum(design["groupsize"]) < mingroupsize):
		raise Exception("sum of groupsizes is less than mintotgroupsize")
	if any(np.sum(design["groupsize"]) > maxgroupsize):
		raise Exception("sum of groupsizes is greater than maxtotgroupsize")
	
	# maxxt and minxt
	if maxxt.size == 1:
		maxxt = ones(size(design["xt"], 1), size(design["xt"], 2)) * maxxt
	if type(maxxt) is list:
		length = max([len(i) for i in maxxt])
		maxxt_ = []
		for i in range(0, len(maxxt)):
			maxxt[i] = maxxt[i].astype(np.float32)
			maxxt[i] = np.pad(maxxt[i], (0, length-len(maxxt[i])), "constant", constant_values=np.nan)
			maxxt_.append(maxxt[i].tolist())
		maxxt = np.array(maxxt_)
	if size(maxxt, 1) == 1 and design["m"] != 1:
		maxxt = np.tile(maxxt.flatten(), design["m"]).reshape(design["m"], maxxt.size)
	if type(maxxt) is not np.ndarray:
		maxxt = np.array(maxxt)
	if size(maxxt, 1) != design["m"]:
		raise Exception("The number of rows in maxxt (" +
						str(size(maxxt, 1)) +
						") is not the same as the number of groups m (" +
						str(design["m"]) +
						")")
	if size(maxxt, 2) == np.max(design["ni"]) and np.max(maxni) > np.max(design["ni"]) and size(design["xt"], 2) == np.max(maxni):
		maxxt_full = design["xt"]
		maxxt_full[:, 0:np.max(design["ni"])] = maxxt
		maxxt = maxxt_full
	if test_mat_size(np.array(design["xt"].shape), np.array(maxxt), "maxxt") == 1:
		maxxt = pd.DataFrame(maxxt,
							 index=["grp_"+str(i) for i in range(1, design["m"]+1)],
							 columns=["obs_"+str(i) for i in range(1, maxxt.shape[1]+1)])

	if minxt.size == 1:
		minxt = ones(size(design["xt"], 1), size(design["xt"], 2)) * minxt
	if type(minxt) is list:
		length = max([len(i) for i in minxt])
		minxt_ = []
		for i in range(0, len(minxt)):
			minxt[i] = minxt[i].astype(np.float32)
			minxt[i] = np.pad(minxt[i], (0, length-len(minxt[i])), "constant", constant_values=np.nan)
			minxt_.append(minxt[i].tolist())
		minxt = np.array(minxt_)
	if size(minxt, 1) == 1 and design["m"] != 1:
		minxt = np.tile(minxt.flatten(), design["m"]).reshape(design["m"], minxt.size)
	if type(minxt) is not np.ndarray:
		minxt = np.array(minxt)
	if size(minxt, 1) != design["m"]:
		raise Exception("The number of rows in minxt (" +
						str(size(minxt, 1)) +
						") is not the same as the number of groups m (" +
						str(design["m"]) +
						")")
	if size(minxt, 2) == np.max(design["ni"]) and np.max(maxni) > np.max(design["ni"]) and size(design["xt"], 2) == np.max(maxni):
		minxt_full = design["xt"]
		minxt_full[:, 0:np.max(design["ni"])] = minxt
		minxt = minxt_full
	if test_mat_size(np.array(design["xt"].shape), np.array(minxt), "minxt") == 1:
		minxt = pd.DataFrame(minxt,
							 index=["grp_"+str(i) for i in range(1, design["m"]+1)],
							 columns=["obs_"+str(i) for i in range(1, minxt.shape[1]+1)])

	# make sure min is less than max
	ret = comp_max_min(maxxt, minxt, called_args)
	maxxt = ret["max_val"]
	minxt = ret["min_val"]

	# check for zeros
	if our_zero is not None:
		minxt = minxt + our_zero * (minxt == 0)
		maxxt = maxxt + our_zero * (maxxt == 0)
		design["xt"] = design["xt"] + our_zero * (design["xt"] == 0)
	
	# check given max and min
	if any(np.greater(minxt, design["xt"])):
		raise Exception("xt is less than minxt")
	if any(np.greater(design["xt"], maxxt)):
		raise Exception("xt is greater than maxxt")
	
	# need to decide on appripriate values of xt and minxt and maxxt if applicable
	if any(maxni>design["ni"]) and any(np.isnan(design["xt"])):
		for grp in range(1, design["m"]+1):
			xt_grp = design["xt"][grp, :]
			maxxt_grp = maxxt[grp, :]
			minxt_grp = minxt[grp, :]
			if any(np.isnan(maxxt_grp)):
				_, idx = np.unique(maxxt_grp[~np.isnan(maxxt_grp)], return_index=True) # remove duplicated and nan values but keep order
				max_vals = [maxxt_grp[~np.isnan(maxxt_grp)][index] for index in sorted(idx)]
				if len(max_vals) == 1:
					maxxt_grp[np.isnan(maxxt_grp)] = max_vals
				else:
					raise Exception("Unable to determine the maxxt values needed for group " +
									str(grp) +
									"\n if ni increases with optimization \nPlease supply them as input.")
			if any(np.isnan(minxt_grp)):
				_, idx = np.unique(minxt_grp[~np.isnan(minxt_grp)], return_index=True) # remove duplicated and nan values but keep order
				min_vals = [minxt_grp[~np.isnan(minxt_grp)][index] for index in sorted(idx)]
				if len(min_vals) == 1:
					minxt_grp[np.isnan(minxt_grp)] = min_vals
				else:
					raise Exception("Unable to determine the minxt values needed for group " +
									str(grp) +
									"\n if ni increases with optimization \nPlease supply them as input.")
			if any(np.isnan(xt_grp)):
				_, idx = np.unique(maxxt_grp[~np.isnan(maxxt_grp)], return_index=True) # remove duplicated and nan values but keep order
				max_vals = [maxxt_grp[~np.isnan(maxxt_grp)][index] for index in sorted(idx)]
				_, idx = np.unique(minxt_grp[~np.isnan(minxt_grp)], return_index=True) # remove duplicated and nan values but keep order
				min_vals = [minxt_grp[~np.isnan(minxt_grp)][index] for index in sorted(idx)]
				one_max = len(max_vals) == 1
				one_min = len(min_vals) == 1
				if one_max and one_min:
					xt_grp[np.isnan(xt_grp)] = np.mean(np.array([max_vals] + [min_vals]))
				if one_max and (not one_min):
					xt_grp[np.isnan(xt_grp)] = max_vals
				if (not one_max) and one_min:
					xt_grp[np.isnan(xt_grp)] = min_vals
				if (not one_max) and (not one_min):
					raise Exception("Unable to determine the initial xt values needed for group " + str(grp) +
                                    "\n if ni increases with optimization \nPlease supply them as input.")
			design["xt"][grp, :] = xt_grp
			maxxt[grp, :] = maxxt_grp
			minxt[grp, :] = minxt_grp
		design_new["xt"] = design["xt"]

	# for a 没写
	if maxa is not None:
		if type(maxa) is list:
			maxa

	# for x
	if x_space is None and design["x"] is not None:
		x_space = cell(size(design["x"]))
		for i in range(1, size(design["x"], 1)+1):
			for j in range(1, size(design["x"], 2)):
				x_space[i, j] = design["x"][i, j]
	if x_space is not None:
		if size(x_space, 1) == 1 and design["m"] != 1:
			x_space = np.tile(x_space.flatten(), design["m"]).reshape(design["m"], x_space.size)
		if test_mat_size(size(design["x"]), x_space, "x_space") == 1:
			x_space = pd.DataFrame(x_space,
								   index=["grp_"+str(i) for i in range(1, design["m"]+1)],
								   columns=design["x"].columns.values.tolist())
		design_space["x_space"] = x_space

		for i in range(1, size(design["x"], 1)+1):
			for j in range(1, size(design["x"], 2)):
				if design["x"][i, j] not in [x_space[i, j]]:
					raise Exception("x value for group " + str(i) + " (column " + str(j) + ") is not in the design space")
	
	# for xt_space






def comp_max_min(max_val, min_val, called_args):
	args = inspect.getargspec(comp_max_min).args
	if any(np.greater(min_val, max_val)):
		 min_val_sup = args[3] in called_args
		 ...
	...