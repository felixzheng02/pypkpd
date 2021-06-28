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
#' @param use_grouped_xt Group sampling times between groups so that each group has the same values (\code{True} or \code{False}).
#' @param grouped_xt Matrix defining the grouping of sample points. Matching integers mean that the points are matched.  
#' Allows for finer control than \code{use_grouped_xt}
#' @param use_grouped_a Group continuous design variables between groups so that each group has the same values (\code{True} or \code{False}).
#' @param grouped_a Matrix defining the grouping of continuous design variables. Matching integers mean that the values are matched.  
#' Allows for finer control than \code{use_grouped_a}.
#' @param use_grouped_x Group discrete design variables between groups so that each group has the same values (\code{True} or \code{False}).
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
						use_grouped_xt=False,
						grouped_xt=None,
						use_grouped_a=False,
						grouped_a=None,
						use_grouped_x=False,
						grouped_x=None,
						our_zero=None):

	called_args = locals()

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
	if size(maxni)[0] == 1 and design["m"] != 1:
		maxni = np.array([maxni] * design["m"]).reshape([design["m"], 1])
	if type(maxni) is not np.ndarray:
		maxni = np.array(maxni)
	if test_mat_size(np.array([design["m"], 1]), maxni, "maxni") == 1:
		maxni = pd.DataFrame(maxni,
							 index=["grp_"+str(i+1) for i in range(0, design["m"])],
							 columns=["n_obs"] * maxni.shape[1])

	# minni
	if size(minni)[0] == 1 and design["m"] != 1:
		minni = np.array([minni] * design["m"]).reshape([design["m"], 1])
	if type(minni) is not np.ndarray:
		minni = np.array(minni)
	if test_mat_size(np.array([design["m"], 1]), minni, "minni") == 1:
		minni = pd.DataFrame(minni,
							 index=["grp_"+str(i+1) for i in range(0, design["m"])],
							 columns=["n_obs"] * minni.shape[1])
	
	# make sure min is smaller than max
	ret = comp_max_min(maxni, minni, called_args)
	maxni = ret["max_val"]
	minni = ret["min_val"]

	# check ni given max and min
	if any(np.array(design["ni"]) < np.array(minni)):
		raise Exception("ni is less than minni")
	if any(np.array(design["ni"]) > np.array(maxni)):
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
	if np.amax(np.array(maxni)) > size(design["xt"])[1]:

		# xt has to increase
		xt_full = np.array(ones(design["m"], int(np.amax(maxni)))) * np.nan
		xt_full[0:design["m"], 0:size(design["xt"])[1]] = design["xt"]
		xt_full = pd.DataFrame(xt_full,
							   index=["grp_"+str(i+1) for i in range(0, design["m"])],
							   columns=["obs_"+str(i+1) for i in range(0, size(xt_full)[1])])
		design["xt"] = xt_full
		design_new["xt"] = design["xt"]

		# model switch has to increase
		model_switch_full = np.array(ones(design["m"], int(np.amax(maxni)))) * np.nan
		model_switch_full[0:design["m"], 0:size(design["model_switch"])[1]] = design["model_switch"]
		model_switch_full = pd.DataFrame(model_switch_full,
							   index=["grp_"+str(i+1) for i in range(0, design["m"])],
							   columns=["obs_"+str(i+1) for i in range(0, size(model_switch_full)[1])])
		design["model_switch"] = model_switch_full
		for i in range(0, design["model_switch"].shape[0]):
			x_tmp = design["model_switch"].iloc[i, :]
			_, idx = np.unique(x_tmp[~np.isnan(x_tmp)], return_index=True) # remove duplicated and nan values but keep order
			if len([x_tmp[~np.isnan(x_tmp)][index] for index in sorted(idx)]) == 1:
				x_tmp[np.isnan(x_tmp)] = [x_tmp[~np.isnan(x_tmp)][index] for index in sorted(idx)]
			else:
				raise Exception("Unable to determine the model_switch values needed for group " + str(i+1)
								+ "\n Please supply them as input.")
			design["model_switch"].iloc[i, :] = x_tmp
		design_new["model_switch"] = design["model_switch"]

	# maxgroupsize
	if size(maxgroupsize)[0] == 1 and design["m"] != 1:
		maxgroupsize = pd.DataFrame(np.array([maxgroupsize] * design["m"]).reshape([design["m"], 1]),
									index=["grp_"+str(i+1) for i in range(0, design["m"])])
	if len(maxgroupsize.shape) != 2:
		maxgroupsize = maxgroupsize[:, np.newaxis]
	if test_mat_size(np.array([design["m"], 1]), np.array(maxgroupsize), "maxgroupsize") == 1:
		maxgroupsize = pd.DataFrame(maxgroupsize,
								 index=["grp_"+str(i) for i in range(1, design["m"]+1)],
								 columns=["n_id"] * maxgroupsize.shape[1])
	
	# mingroupsize
	if size(mingroupsize)[0] == 1 and design["m"] != 1:
		mingroupsize = pd.DataFrame(np.array([mingroupsize] * design["m"]).reshape([design["m"], 1]),
									index=["grp_"+str(i+1) for i in range(0, design["m"])])
	if len(mingroupsize.shape) != 2:
		mingroupsize = mingroupsize[:, np.newaxis]
	if test_mat_size(np.array([design["m"], 1]), np.array(mingroupsize), "mingroupsize") == 1:
		mingroupsize = pd.DataFrame(mingroupsize,
								 index=["grp_"+str(i) for i in range(1, design["m"]+1)],
								 columns=["n_id"] * mingroupsize.shape[1])

	# make sure min is less than max
	ret = comp_max_min(maxgroupsize, mingroupsize, called_args)
	maxgroupsize = ret["max_val"]
	mingroupsize = ret["min_val"]

	# check given max and min
	if any(np.array(design["groupsize"]) < np.array(mingroupsize)):
		raise Exception("groupsize is less than mingroupsize")
	if any(np.array(design["groupsize"]) > np.array(maxgroupsize)):
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
	if any(np.array(np.sum(design["groupsize"])) < np.array(mintotgroupsize)):
		raise Exception("sum of groupsizes is less than mintotgroupsize")
	if any(np.array(np.sum(design["groupsize"])) > np.array(maxtotgroupsize)):
		raise Exception("sum of groupsizes is greater than maxtotgroupsize")
	
	# maxxt and minxt
	if type(maxxt) is int:
		maxxt = np.array(ones(size(design["xt"])[0], size(design["xt"])[1])) * maxxt
	elif maxxt.size == 1:
		maxxt = np.array(ones(size(design["xt"])[0], size(design["xt"])[1])) * maxxt
	if type(maxxt) is list:
		length = max([len(i) for i in maxxt])
		maxxt_ = []
		for i in range(0, len(maxxt)):
			maxxt[i] = maxxt[i].astype(np.float32)
			maxxt[i] = np.pad(maxxt[i], (0, length-len(maxxt[i])), "constant", constant_values=np.nan)
			maxxt_.append(maxxt[i].tolist())
		maxxt = np.array(maxxt_)
	if size(maxxt)[0] == 1 and design["m"] != 1:
		maxxt = np.tile(maxxt.flatten(), design["m"]).reshape(design["m"], maxxt.size)
	if type(maxxt) is not np.ndarray:
		maxxt = np.array(maxxt)
	if size(maxxt)[0] != design["m"]:
		raise Exception("The number of rows in maxxt (" +
						str(size(maxxt)[0]) +
						") is not the same as the number of groups m (" +
						str(design["m"]) +
						")")
	if size(maxxt)[1] == int(np.max(design["ni"])) and int(np.max(maxni)) > int(np.max(design["ni"])) and size(design["xt"])[1] == int(np.max(maxni)):
		maxxt_full = design["xt"]
		maxxt_full[:, 0:np.max(design["ni"])] = maxxt
		maxxt = maxxt_full
	if test_mat_size(np.array(design["xt"].shape), np.array(maxxt), "maxxt") == 1:
		maxxt = pd.DataFrame(maxxt,
							 index=["grp_"+str(i+1) for i in range(0, design["m"])],
							 columns=["obs_"+str(i+1) for i in range(0, maxxt.shape[1])])

	if type(minxt) is int:
		minxt = np.array(ones(size(design["xt"])[0], size(design["xt"])[1])) * minxt
	elif minxt.size == 1:
		minxt = np.array(ones(size(design["xt"])[0], size(design["xt"])[1])) * minxt
	if type(minxt) is list:
		length = max([len(i) for i in minxt])
		minxt_ = []
		for i in range(0, len(minxt)):
			minxt[i] = minxt[i].astype(np.float32)
			minxt[i] = np.pad(minxt[i], (0, length-len(minxt[i])), "constant", constant_values=np.nan)
			minxt_.append(minxt[i].tolist())
		minxt = np.array(minxt_)
	if size(minxt)[0] == 1 and design["m"] != 1:
		minxt = np.tile(minxt.flatten(), design["m"]).reshape(design["m"], minxt.size)
	if type(minxt) is not np.ndarray:
		minxt = np.array(minxt)
	if size(minxt)[0] != design["m"]:
		raise Exception("The number of rows in minxt (" +
						str(size(minxt)[0]) +
						") is not the same as the number of groups m (" +
						str(design["m"]) +
						")")
	if size(minxt)[1] == int(np.max(design["ni"])) and int(np.max(maxni)) > int(np.max(design["ni"])) and size(design["xt"])[1] == int(np.max(maxni)):
		minxt_full = design["xt"]
		minxt_full[:, 0:np.max(design["ni"])] = minxt
		minxt = minxt_full
	if test_mat_size(np.array(design["xt"].shape), np.array(minxt), "minxt") == 1:
		minxt = pd.DataFrame(minxt,
							 index=["grp_"+str(i+1) for i in range(0, design["m"])],
							 columns=["obs_"+str(i+1) for i in range(0, minxt.shape[1])])

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
	if np.greater(np.array(minxt), np.array(design["xt"])).any():
		raise Exception("xt is less than minxt")
	if np.greater(np.array(design["xt"]), np.array(maxxt)).any():
		raise Exception("xt is greater than maxxt")
	
	# need to decide on appripriate values of xt and minxt and maxxt if applicable
	if any(np.array(maxni) > np.array(design["ni"])) and any(np.isnan(design["xt"])):
		for grp in range(0, design["m"]):
			xt_grp = design["xt"].iloc[grp, :]
			maxxt_grp = pd.DataFrame(maxxt).iloc[grp, :]
			minxt_grp = pd.DataFrame(minxt).iloc[grp, :]
			if any(np.isnan(maxxt_grp)):
				_, idx = np.unique(maxxt_grp[~np.isnan(maxxt_grp)], return_index=True) # remove duplicated and nan values but keep order
				max_vals = [maxxt_grp[~np.isnan(maxxt_grp)][index] for index in sorted(idx)]
				if len(max_vals) == 1:
					maxxt_grp[np.isnan(maxxt_grp)] = max_vals
				else:
					raise Exception("Unable to determine the maxxt values needed for group " +
									str(grp+1) +
									"\n if ni increases with optimization \nPlease supply them as input.")
			if any(np.isnan(minxt_grp)):
				_, idx = np.unique(minxt_grp[~np.isnan(minxt_grp)], return_index=True) # remove duplicated and nan values but keep order
				min_vals = [minxt_grp[~np.isnan(minxt_grp)][index] for index in sorted(idx)]
				if len(min_vals) == 1:
					minxt_grp[np.isnan(minxt_grp)] = min_vals
				else:
					raise Exception("Unable to determine the minxt values needed for group " +
									str(grp+1) +
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
					raise Exception("Unable to determine the initial xt values needed for group " + str(grp+1) +
                                    "\n if ni increases with optimization \nPlease supply them as input.")
			design["xt"].iloc[grp, :] = xt_grp
			pd.DataFrame(maxxt).iloc[grp, :] = maxxt_grp
			pd.DataFrame(minxt).iloc[grp, :] = minxt_grp
		design_new["xt"] = design["xt"]

	# for a
	if maxa is not None:
		if type(maxa) is list:
			maxa = pd.DataFrame(maxa)
		if size(maxa)[0] == 1 and design["m"] != 1:
			maxa = np.tile(maxa.flatten(), design["m"]).reshape(design["m"], maxa.size)
		if type(maxa) is not np.ndarray and type(maxa) is not pd.DataFrame:
			maxa = pd.DataFrame(maxa)
		if size(maxa)[0] != design["m"]:
			raise Exception("The number of rows in maxa (" +
							str(size(maxa)[0]) +
							") is not the same as the number of groups m (" +
							str(design["m"]) + ")")
		maxa.set_axis(["grp_" + str(i+1) for i in range(0, design["m"])])
		if maxa.columns.values.tolist() == []:
			maxa.set_axis(design["a"].columns.values.tolist())
		design_space["maxa"] = maxa
	
	if mina is not None:
		if type(mina) is list:
			mina = pd.DataFrame(mina)
		if size(mina)[0] == 1 and design["m"] != 1:
			mina = np.tile(mina.flatten(), design["m"]).reshape(design["m"], mina.size)
		if type(mina) is not np.ndarray and type(mina) is not pd.DataFrame:
			mina = pd.DataFrame(mina)
		if size(mina)[0] != design["m"]:
			raise Exception("The number of rows in maxa (" +
							str(size(mina)[0]) +
							") is not the same as the number of groups m (" +
							str(design["m"]) + ")")
		mina.set_axis(["grp_" + str(i+1) for i in range(0, design["m"])])
		if mina.columns.values.tolist() == []:
			mina.set_axis(design["a"].columns.values.tolist())
		design_space["mina"] = mina

	# make sure max is min smaller than max
	if mina is not None and maxa is not None:
		ret = comp_max_min(maxa, mina, called_args)
		maxa = ret["max_val"]
		mina = ret["min_val"]

	# check ni given max and min
	if mina is not None and maxa is not None and design["a"] is not None:
		if (np.array(design["a"]) < np.array(mina)).any():
			raise Exception("a is less than mina")
		if (np.array(design["a"]) > np.array(maxa)).any():
			raise Exception("a is greater than maxa")

	# for x
	if x_space is None and "x" in design.keys():
		x_space = cell(size(design["x"]))
		for i in range(0, size(design["x"])[0]):
			for j in range(0, size(design["x"])[1]):
				x_space[i, j] = design["x"][i, j]
	if x_space is not None:
		if size(x_space)[0] == 1 and design["m"] != 1:
			x_space = np.tile(x_space.flatten(), design["m"]).reshape(design["m"], x_space.size)
		if test_mat_size(np.array(size(design["x"])), x_space, "x_space") == 1:
			x_space = pd.DataFrame(x_space,
								   index=["grp_"+str(i+1) for i in range(0, design["m"])],
								   columns=design["x"].columns.values.tolist())
		design_space["x_space"] = x_space

		for i in range(0, size(design["x"])[0]):
			for j in range(0, size(design["x"])[1]):
				if design["x"][i, j] not in [x_space[i, j]]:
					raise Exception("x value for group " + str(i+1) + " (column " + str(j+1) + ") is not in the design space")
	
	# for xt_space
	if xt_space is not None:
		if type(xt_space) is list:  # then it is a list with no dimensions, need to convert to a cell
			nspace = len(xt_space)
			nrow_xt = design["xt"].shape[0]
			ncol_xt = design["xt"].shape[1]
			if nspace == 1: # all time points in all groups have the same space
				xt_space_tmp = design["xt"]
				xt_space = cell(size(design["xt"]))
				xt_space = xt_space_tmp
			elif nspace == ncol_xt: # we assume that all groups have the same space
				xt_space_tmp = xt_space
				xt_space = cell(size(design["xt"]))
				for i in range(0, nrow_xt):
					xt_space[i, :] = xt_space_tmp
			elif nspace == (ncol_xt * nrow_xt): # we assume that spaces are entered in row major form
				xt_space_tmp = xt_space
				xt_space = np.array(xt_space_tmp).reshape([int((xt_space_tmp.size)/nrow_xt), nrow_xt])
		else: # then it is a vector, assume the vector is the same for all xt's
			tmp_lst = xt_space.tolist()
			xt_space = cell(size(design["xt"]))
			xt_space = tmp_lst
		if size(xt_space)[0] == 1 and design["m"] != 1:
			xt_space = np.tile(xt_space.flatten(), design["m"]).reshape(design["m"], xt_space.size)
		if size(xt_space)[1] == 1 and size(design["xt"])[2] != 1:
			xt_space = np.transpose(np.tile(xt_space.flatten(), size(design["xt"])[1]).reshape(design["m"], size(design["xt"])[1]))
		if test_mat_size(np.array(size(design["xt"])), np.array(xt_space), "xt_space") == 1:
			xt_space = pd.DataFrame(xt_space,
									index=["grp_"+str(i+1) for i in range(0, design["m"])],
									columns=design["xt"].columns.values.tolist())
		
		for i in range(0, design["xt"].shape[0]):
			for j in range(0, design["xt"].shape[1]):
				if design["xt"][i, j] not in [xt_space[i, j]] and design["xt"][i, j] != np.nan:
					raise Exception("xt value for group " + str(i+1) + " (column " + str(j+1) + ") is not in the design space")

	# for a_space
	if a_space is not None:
		if type(a_space) is not np.ndarray and type(a_space) is list:
			a_space = np.array(a_space * design["m"]).reshape([design["m"], len(a_space)])
		else:
			tmp_lst = pd.DataFrame(np.array(a_space).reshape([1, a_space.size]),
								   # 没写 dimnames = list(None,names(x))
								   )
			mat = pd.DataFrame([])
			for jj in range(0, tmp_lst.size):
				tmp = tmp_lst[0, jj]
				if any(type(i) == str for i in tmp.columns.values.tolist()) and any(type(i) == str for i in design["a"].columns.values.tolist()):
					tmp = tmp[design["a"].columns.values.tolist()]
					mat = mat.append(tmp)
			a_space = mat
		if size(a_space)[0] == 1 and design["m"] != 1:
			a_space = np.tile(np.array(a_space).flatten(), design["m"]).reshape(design["m"], a_space.size)
		if size(a_space)[1] == 1 and size(design["a"])[1] == 1:
			a_space = np.transpose(np.tile(a_space.flatten(), size(design["a"])[1]).reshape(design["m"], size(design["a"])[1]))
		if test_mat_size(np.array(size(design["a"])), np.array(a_space), "a_space") == 1:
			if type(a_space) is not np.ndarray and type(a_space) is not pd.DataFrame and all(np.array(size(a_space)) == 1):
				a_space = np.array(a_space)
				a_space = pd.DataFrame(a_space,
									   index=["grp_" + str(i+1) for i in range(0, design["m"])],
									   columns=design["a"].columns.values.tolist())
		for i in range(0, design["a"].shape[0]):
			for j in range(0, design["a"].shape[1]):
				if design["a"][i, j] not in [a_space[i, j]] and design["a"] is not None:
					raise Exception("a value for group " + str(i+1) + " (column " + str(j+1) + ") is not in the design space")

	# for grouped_xt
	if grouped_xt is None:
		grouped_xt = design["xt"] * np.nan
		val = 1
		for i in range(0, design["xt"].shape[0]):
			if use_grouped_xt:
				val = 1
			for j in range(0, design["xt"].shape[1]):
				if design["xt"].iloc[i, j] is not None:
					grouped_xt.iloc[i, j] = val
					val += 1

	if grouped_xt.size == 1:
		grouped_xt = np.array(ones(size(design["xt"])[0], size(design["xt"])[1])) * grouped_xt
		use_grouped_xt = True
	if type(grouped_xt) is list:
		length = max([len(i) for i in grouped_xt])
		grouped_xt_ = []
		for i in range(0, len(grouped_xt)):
			grouped_xt[i] = grouped_xt[i].astype(np.float32)
			grouped_xt[i] = np.pad(grouped_xt[i], (0, length-len(grouped_xt[i])), "constant", constant_values=np.nan)
			grouped_xt_.append(grouped_xt[i].tolist())
		grouped_xt = np.array(grouped_xt_)
	if size(grouped_xt)[0] == 1 and design["m"] != 1:
		grouped_xt = np.tile(grouped_xt.flatten(), design["m"]).reshape(design["m"], grouped_xt.size)
		use_grouped_xt = True
	if type(grouped_xt) is not np.ndarray and type(grouped_xt) is not pd.DataFrame:
		grouped_xt = np.array(grouped_xt)
	if size(grouped_xt)[1] == size(grouped_xt)[1] and np.max(np.array(maxni)) > np.max(np.array(design["ni"])) and size(design["xt"])[1] == np.max(np.array(maxni)):
		grouped_xt_full = design["xt"] * np.nan
		grouped_xt_full.iloc[:, 0:int(np.max(design["ni"]))] = pd.DataFrame(grouped_xt.astype(np.float32))
		grouped_xt = grouped_xt_full
	if test_mat_size(np.array(size(design["xt"])), np.array(grouped_xt), "grouped_xt") == 1:
		grouped_xt = pd.DataFrame(grouped_xt,
								  index=["grp_" + str(i+1) for i in range(0, design["m"])],
								  columns=["obs_" + str(i+1) for i in range(0, grouped_xt.shape[1])])

	# get values in the NA region if possible
	if any(maxni > design["ni"]) and any(np.isnan(grouped_xt)):
		for grp in range(0, design["m"]):
			grouped_xt_grp = pd.DataFrame(grouped_xt).iloc[grp, :]
			if any(np.isnan(grouped_xt_grp)):
				_, idx = np.unique(grouped_xt_grp[~np.isnan(grouped_xt_grp)], return_index=True) # remove duplicated and nan values but keep order
				vals = [grouped_xt_grp[~np.isnan(grouped_xt_grp)][index] for index in sorted(idx)]
				if len(vals) == 1:
					grouped_xt_grp[np.isnan(grouped_xt_grp)] = vals
				else:
					raise Exception("Unable to determine the grouped_xt values needed for group " +
									str(grp+1) +
									"\n if ni increases with optimization \nPlease supply them as input.")
			pd.DataFrame(grouped_xt).iloc[grp, :] = grouped_xt_grp

	_, idx = np.unique(np.array(grouped_xt)[~np.isnan(np.array(design["xt"]))], return_index=True) # remove duplicated and nan values but keep order
	for i in [np.array(grouped_xt)[~np.isnan(np.array(design["xt"]))][index] for index in sorted(idx)]:
		_, idx = np.unique(np.array(design["xt"])[np.logical_and(np.array(grouped_xt) == i, ~np.isnan(np.array(design["xt"])))], return_index=True) # remove duplicated and nan values but keep order
		if len([np.array(design["xt"])[np.logical_and(np.array(grouped_xt) == i, ~np.isnan(np.array(design["xt"])))][index] for index in sorted(idx)]) != 1:
			raise Exception("xt values grouped with value %g from grouped_xt do not have the same initial values.\n'" % i)
		_, idx = np.unique(np.array(maxxt)[np.logical_and(np.array(grouped_xt) == i, ~np.isnan(np.array(design["xt"])))], return_index=True) # remove duplicated and nan values but keep order
		if len([np.array(maxxt)[np.logical_and(np.array(grouped_xt) == i, ~np.isnan(np.array(design["xt"])))][index] for index in sorted(idx)]) != 1:
			raise Exception("xt values grouped with value %g from grouped_xt do not have the same maximum allowed values (maxxt).\n" % i)
		_, idx = np.unique(np.array(minxt)[np.logical_and(np.array(grouped_xt) == i, ~np.isnan(np.array(design["xt"])))], return_index=True) # remove duplicated and nan values but keep order
		if len([np.array(minxt)[np.logical_and(np.array(grouped_xt) == i, ~np.isnan(np.array(design["xt"])))][index] for index in sorted(idx)]) != 1:
			raise Exception("xt values grouped with value %g from grouped_xt do not have the same maximum allowed values (minxt).\n" % i)
		if xt_space is None:
			grouped_cells_xt = None
		else:
			grouped_cells_xt = np.array(xt_space)[np.logical_and(np.array(grouped_xt) == i, ~np.isnan(np.array(design["xt"])))]
			for j in range(0, grouped_cells_xt.size):
				for k in range(0, grouped_cells_xt.size):
					if (any(np.array(size(grouped_cells_xt[j])) != np.array(size(grouped_cells_xt[k]))) or
						(grouped_cells_xt[j] != grouped_cells_xt[k])):
						raise Exception("xt values grouped with value %g from grouped_xt do not have the same allowed discrete values (xt_space).\n" % i)

	_, idx = np.unique(np.array(grouped_xt)[~np.isnan(np.array(design["xt"]))], return_index=True) # remove duplicated and nan values but keep order
	for i in range(0, int(np.max([np.array(grouped_xt)[~np.isnan(np.array(design["xt"]))][index] for index in sorted(idx)]))):
		_, idx = np.unique(np.array(design["xt"])[np.logical_and(np.array(grouped_xt) == i+1, ~np.isnan(np.array(design["xt"])))], return_index=True) # remove duplicated and nan values but keep order
		if len([np.array(design["xt"])[np.logical_and(np.array(grouped_xt) == i+1, ~np.isnan(np.array(design["xt"])))][index] for index in sorted(idx)]) == 0:
			raise Exception("grouped_xt must be sequential and cannot have missing values.\nNo xt values were grouped with value %g in grouped_xt.\n" % i)

	# for grouped_a
	if "a" in design.keys():
		if design["a"] is not None:
			if grouped_a is None:
				grouped_a = np.array(design["a"]) * np.nan
				val = 1
				for i in range(0, size(design["a"])[0]):
					if use_grouped_a:
						val = 1
					for j in range(0, size(design["a"])[1]):
						if ~np.isnan(np.array(design["a"])[i, j]):
							grouped_a[i, j] = val
							val += 1

			if grouped_a.size == 1:
				grouped_a = np.array(ones(size(design["a"][0]), size(design["a"][1]))) * grouped_a

			if type(grouped_a) == list:
				length = max([len(i) for i in grouped_a])
				grouped_a_ = []
				for i in range(0, len(grouped_a)):
					grouped_a[i] = grouped_a[i].astype(np.float32)
					grouped_a[i] = np.pad(grouped_a[i], (0,length-len(grouped_a[i])),"constant",constant_values=np.nan)
					grouped_a_.append(grouped_a[i].tolist())
				grouped_a = np.array(grouped_a_)
				use_grouped_a = True

			if type(grouped_a) is not np.ndarray and type(grouped_a) is not pd.DataFrame:
				grouped_a = np.array(grouped_a)

			if test_mat_size(np.array(size(design["a"])), grouped_a, "grouped_a") == 1:
				grouped_a = pd.DataFrame(grouped_a,
										index=["grp_" + str(i+1) for i in range(0, design["m"])])
				if grouped_a.columns.values.tolist() == []:
					grouped_a.set_axis(design["a"].columns.values.tolist(), axis="column")

			_, idx = np.unique(np.array(grouped_a)[~np.isnan(np.array(design["a"]))], return_index=True) # remove duplicated and nan values but keep order
			for i in [np.array(grouped_a)[~np.isnan(np.array(design["a"]))][index] for index in sorted(idx)]:
				_, idx = np.unique(np.array(design["a"])[np.logical_and(np.array(grouped_a) == i, ~np.isnan(np.array(design["a"])))], return_index=True) # remove duplicated and nan values but keep order
				if len([np.array(design["a"])[np.logical_and(np.array(grouped_a) == i, ~np.isnan(np.array(design["a"])))][index] for index in sorted(idx)]) != 1:
					raise Exception("a values grouped with value %g from grouped_a do not have the same initial values.\n'" % i)
				_, idx = np.unique(np.array(maxa)[np.logical_and(np.array(grouped_a) == i, ~np.isnan(np.array(design["a"])))], return_index=True) # remove duplicated and nan values but keep order
				if len([np.array(maxa)[np.logical_and(np.array(grouped_a) == i, ~np.isnan(np.array(design["a"])))][index] for index in sorted(idx)]) != 1:
					raise Exception("a values grouped with value %g from grouped_a do not have the same maximum allowed values (maxa).\n" % i)
				_, idx = np.unique(np.array(mina)[np.logical_and(np.array(grouped_a) == i, ~np.isnan(np.array(design["a"])))], return_index=True) # remove duplicated and nan values but keep order
				if len([np.array(mina)[np.logical_and(np.array(grouped_a) == i, ~np.isnan(np.array(design["a"])))][index] for index in sorted(idx)]) != 1:
					raise Exception("a values grouped with value %g from grouped_a do not have the same maximum allowed values (mina).\n" % i)
				if a_space is None:
					grouped_cells_a = None
				else:
					grouped_cells_a = a_space[np.logical_and(np.array(grouped_a) == i, np.isnan(np.array(np.array(design["a"]))))]
					for j in range(0, grouped_cells_a.size):
						for k in range(0, grouped_cells_a.size):
							if (any(np.array(size(grouped_cells_a[j])) != np.array(size(grouped_cells_a[k]))) or
								(grouped_cells_a[j] != grouped_cells_a[k])): # 没写any
								raise Exception("a values grouped with value %g from grouped_a do not have the same allowed discrete values (a_space).\n" % i)

			_, idx = np.unique(np.array(grouped_a)[~np.isnan(np.array(design["a"]))], return_index=True) # remove duplicated and nan values but keep order
			for i in range(0, int(np.max([np.array(grouped_a)[~np.isnan(np.array(design["a"]))][index] for index in sorted(idx)]))):
				_, idx = np.unique(np.array(design["a"])[np.logical_and(np.array(grouped_a) == i+1, ~np.isnan(np.array(design["a"])))], return_index=True) # remove duplicated and nan values but keep order
				if len([np.array(design["a"])[np.logical_and(np.array(grouped_a) == i+1, ~np.isnan(np.array(design["a"])))][index] for index in sorted(idx)]) == 0:
					raise Exception("grouped_a must be sequential and cannot have missing values.\nNo a values were grouped with value %g in grouped_a.\n" % i)

			design_space["grouped_a"] = grouped_a
			design_space["use_grouped_a"] = use_grouped_a

	# for grouped_x
	if "x" in design.keys():
		if design["x"] is not None:
			if grouped_x is None:
				grouped_x = design["x"] * np.nan
				val = 1
				for i in range(0, size(design["x"])[0]):
					if use_grouped_x:
						val = 1
					for j in range(0, size(design["x"])[1]):
						if ~np.isnan(design["x"][i, j]):
							grouped_x[i, j] = val
							val += 1

			if grouped_x.size == 1:
				grouped_x = np.array(ones(size(design["x"][0]), size(design["x"][1]))) * grouped_x

			if type(grouped_x) == list:
				length = max([len(i) for i in grouped_x])
				grouped_x_ = []
				for i in range(0, len(grouped_x)):
					grouped_x[i] = grouped_x[i].astype(np.float32)
					grouped_x[i] = np.pad(grouped_x[i], (0,length-len(grouped_x[i])),"constant",constant_values=np.nan)
					grouped_x_.append(grouped_x[i].tolist())
				grouped_x = np.array(grouped_x_)
				use_grouped_x = True

			if type(grouped_x) is not np.ndarray and type(grouped_x) is not pd.DataFrame:
				grouped_x = np.array(grouped_x)

			if test_mat_size(np.array(size(design["x"])), grouped_x, "grouped_x") == 1:
				grouped_x = pd.DataFrame(grouped_x,
										index=["grp_" + str(i+1) for i in range(0, design["m"])])
				if grouped_x.columns.values.tolist() == []:
					grouped_x.set_axis(design["x"].columns.values.tolist(), axis="column")

			_, idx = np.unique(grouped_x[~np.isnan(design["x"])], return_index=True) # remove duplicated and nan values but keep order
			for i in [grouped_x[~np.isnan(design["x"])][index] for index in sorted(idx)]:
				_, idx = np.unique(design["x"][grouped_x == i and ~np.isnan(design["x"])], return_index=True) # remove duplicated and nan values but keep order
				if len([design["x"][grouped_x == i and ~np.isnan(design["x"])][index] for index in sorted(idx)]) != 1:
					raise Exception("x values grouped with value %g from grouped_x do not have the same initial values.\n'" % i)
				grouped_cells = x_space[grouped_x == i and np.isnan(design["x"])]
				for j in range(0, grouped_cells.size):
					for k in range(0, grouped_cells.size):
						if (any(np.array(size(grouped_cells[j])) != np.array(size(grouped_cells[k]))) or
							any(grouped_cells[j] != grouped_cells[k])):
							raise Exception("x values grouped with value %g from grouped_x do not have the same allowed discrete values (grouped_x).\n" % i)

			_, idx = np.unique(grouped_x[~np.isnan(design["x"])], return_index=True) # remove duplicated and nan values but keep order
			for i in range(0, np.max([grouped_x[~np.isnan(design["x"])][index] for index in sorted(idx)])):
				_, idx = np.unique(design["x"][grouped_x == i+1 and ~np.isnan(design["x"])], return_index=True) # remove duplicated and nan values but keep order
				if len([design["x"][grouped_x == i+1 and ~np.isnan(design["x"])][index] for index in sorted(idx)]) == 0:
					raise Exception("grouped_x must be sequential and cannot have missing values.\nNo x values were grouped with value %g in grouped_x.\n" % i)

			design_space["grouped_x"] = grouped_x
			design_space["use_grouped_x"] = use_grouped_x

	design_space["maxni"] = maxni
	design_space["minni"] = minni

	design_space["maxtotni"] = maxtotni
	design_space["mintotni"] = mintotni

	design_space["maxgroupsize"] = maxgroupsize
	design_space["mingroupsize"] = mingroupsize

	design_space["maxtotgroupsize"] = maxtotgroupsize
	design_space["mintotgroupsize"] = mintotgroupsize

	design_space["maxxt"] = maxxt
	design_space["minxt"] = minxt	
	design_space["xt_space"] = xt_space

	design_space["grouped_xt"] = grouped_xt
	design_space["use_grouped_xt"] = use_grouped_xt

	design_space["a_space"] = a_space

	# update max and min of a and xt if imputed and discrete
	if maxa_imputed and a_space is not None:
		for i in range(0, a_space.shape[0]):
			for j in range(0, a_space.shape[1]):
				maxa[i, j] = np.max(a_space[i, j][0])
	if mina_imputed and a_space is not None:
		for i in range(0, a_space.shape[0]):
			for j in range(0, a_space.shape[1]):
				mina[i, j] = np.min(a_space[i, j][0])
	design_space["maxa"] = maxa
	design_space["mina"] = mina

	if maxxt_imputed and xt_space is not None:
		for i in range(0, xt_space.shape[0]):
			for j in range(0, xt_space.shape[1]):
				maxxt[i, j] = np.max(xt_space[i, j][0])
	if minxt_imputed and xt_space is not None:
		for i in range(0, xt_space.shape[0]):
			for j in range(0, xt_space.shape[1]):
				minxt[i, j] = np.min(xt_space[i, j][0])
	design_space["maxxt"] = maxxt
	design_space["minxt"] = minxt

	return {"design": design_new, "design_space": design_space}

def comp_max_min(max_val, min_val, called_args):
	args = list(locals().values())
	if any(np.greater(min_val, max_val)):
		min_val_sup = str(args[1]) in list(called_args.keys())
		max_val_sup = str(args[0]) in list(called_args.keys())
		if min_val_sup and max_val_sup:
			raise Exception("Some value of " + str(args[0]) + " is smaller than " + args[1])
		if min_val_sup and ~max_val_sup:
			max_val = np.maximum(max_val, min_val)
		if ~min_val_sup and max_val_sup:
			min_val = np.minimum(max_val, min_val)
	return {"max_val": max_val, "min_val": min_val}