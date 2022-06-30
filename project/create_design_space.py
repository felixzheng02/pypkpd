"""
## Create design variables and a design space for a full description of an optimization problem.
## 
## \code{create_design_space} takes an initial design and arguments for a design space and 
## creates a design and design space for design optimization.
## Checks the sizes of supplied design space variables and 
## changes them to sizes that  make sense if there are inconsistencies.  
## Function arguments can use shorthand notation (single values, vectors, lists of vectors and 
## list of list) or matricies.
## Returns a list of matricies compatible with PopED.
## 
## If a value (or a vector or a list of values) is supplied that corresponds to only one group and the design has
## multiple groups then all groups will have the same value(s). If a Matrix is expected then a list of lists can be supplied 
## instead, each list corresponding to a group.   
## 
## @param design  The output from a call to \code{\link{create_design}}.
## @param maxni Vector defining the maximum number of samples per group. 
## @param minni Vector defining the minimum number of samples per group. 
## @param maxtotni Number defining the maximum number of samples allowed in the experiment. 
## @param mintotni Number defining the minimum number of samples allowed in the experiment.  
## @param maxgroupsize Vector defining the maximum size of the different groups (maximum number of individuals in each group)
## @param mingroupsize Vector defining the minimum size of the different groups (minimum num individuals in each group) 
## @param maxtotgroupsize The total maximal groupsize over all groups
## @param mintotgroupsize The total minimal groupsize over all groups
## @param maxxt Matrix or single value defining the maximum value for each xt sample.  If a single value is 
## supplied then all xt values are given the same maximum value.
## @param minxt Matrix or single value defining the minimum value for each xt sample.  If a single value is 
## supplied then all xt values are given the same minimum value
## @param x_space Cell array \code{\link{cell}} defining the discrete variables for each x value.
## @param xt_space Cell array \code{\link{cell}} defining the discrete variables allowed for each xt value.
##   Can also be a vector of values \code{c(1:10)} (same values allowed for all xt), or a list of lists 
##  \code{list(1:10, 2:23, 4:6)} (one for each value in xt in row major order or just for one row in xt, 
##  and all other rows will be duplicated).
## @param a_space Cell array \code{\link{cell}} defining the discrete variables allowed for each a value.
##   Can also be a list of values \code{list(1:10)} (same values allowed for all a), or a list of lists 
##  \code{list(1:10, 2:23, 4:6)} (one for each value in a).
## @param maxa Vector defining the maximum value for each covariate. IF a single value is supplied then
##  all a values are given the same maximum value
## @param mina Vector defining the minimum value for each covariate. IF a single value is supplied then
##  all a values are given the same minimum value
## @param use_grouped_xt Group sampling times between groups so that each group has the same values (\code{True} or \code{False}).
## @param grouped_xt Matrix defining the grouping of sample points. Matching integers mean that the points are matched.  
## Allows for finer control than \code{use_grouped_xt}
## @param use_grouped_a Group continuous design variables between groups so that each group has the same values (\code{True} or \code{False}).
## @param grouped_a Matrix defining the grouping of continuous design variables. Matching integers mean that the values are matched.  
## Allows for finer control than \code{use_grouped_a}.
## @param use_grouped_x Group discrete design variables between groups so that each group has the same values (\code{True} or \code{False}).
## @param grouped_x Matrix defining the grouping of discrete design variables. Matching integers mean that the values are matched.  
## Allows for finer control than \code{use_grouped_x}.
## @param our_zero Value to interpret as zero in design.
## 
## 
## @family poped_input
##  
## @export
## 
# @importFrom dplyr rbind_all
Author: Caiya Zhang, Yuchen Zheng
"""


import path
import numpy as np
import copy
from project.test_mat_size import test_mat_size
from project.ones import ones
from project.cell import cell
from matpy.matrix import Matrix
from project.util import default_if_none


def create_design_space(design_,
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

	design = copy.deepcopy(design_)
	called_args = locals()

	# assign defaults if not supplied
	maxni = default_if_none(maxni, design["ni"])
	if type(maxni) is not Matrix:
		maxni = Matrix(maxni)
	minni = default_if_none(minni, design["ni"])
	if type(minni) is not Matrix:
		minni = Matrix(minni)
	maxgroupsize = default_if_none(maxgroupsize, design["groupsize"])
	if type(maxgroupsize) is not Matrix:
		maxgroupsize = Matrix(maxgroupsize)
	mingroupsize = default_if_none(mingroupsize, design["groupsize"])
	if type(mingroupsize) is not Matrix:
		mingroupsize = Matrix(mingroupsize)
	
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

	m: int = design["m"].get_value()
	ni: Matrix = design["ni"]
	xt: Matrix = design["xt"]
	model_switch: Matrix = design["model_switch"]
	groupsize: Matrix = design["groupsize"]
	a: Matrix = design["a"]
	if "x" in design.keys():
		x: Matrix = design["x"]
	else:
		x = None


	# maxni
	if maxni.get_shape()[0] == 1 and m != 1:
		maxni.repeat([1, m], [m, 1])
	if test_mat_size([m, 1], maxni.get_shape(), "maxni") == 1:
		maxni.set_axisnam(["grp_"+str(i+1) for i in range(0, m)], 
							["n_obs"] * maxni.get_shape()[1])				  
	# minni
	if minni.get_shape()[0] == 1 and m != 1:
		minni.repeat([1, m], [m, 1])
	if test_mat_size([m, 1], minni.get_shape(), "minni") == 1:
		minni.set_axisnam(["grp_"+str(i+1) for i in range(0, m)], 
							["n_obs"] * minni.get_shape[1])
	
	# make sure min is smaller than max
	ret = comp_max_min(maxni, minni, called_args)
	maxni = Matrix(ret["max_val"])
	minni = Matrix(ret["min_val"])

	# check ni given max and min
	if any(ni.get_data() < minni.get_data()):
		raise Exception("ni is less than minni")
	if any(ni.get_data() > maxni.get_data()):
		raise Exception("ni is greater than maxni")
	

	# maxtotni and mintotni
	if maxtotni is None:
		maxtotni = Matrix(np.sum(maxni.get_data()))
	if mintotni is None:
		mintotni = Matrix(np.sum(minni.get_data()))
	test_mat_size([1, 1], maxtotni.get_shape(), "maxtotni")
	test_mat_size([1, 1], mintotni.get_shape(), "mintotni")
	ret = comp_max_min(maxtotni.get_data(), mintotni.get_data(), called_args)
	maxtotni = Matrix(ret["max_val"])
	mintotni = Matrix(ret["min_val"])
	if any(ni.get_data() < mintotni.get_data()):
		raise Exception("sum of ni is less than mintotni")
	if any(ni.get_data() > maxtotni.get_data()):
		raise Exception("sum of ni is greater than maxtotni")


	# update xt and model_switch given maxni
	if np.amax(maxni.get_data()) > xt.get_shape()[1]:

		# xt has to increase
		xt_full = Matrix(np.ones(m, int(np.amax(maxni.get_data()))) * np.nan)
		xt_full.get_data()[0:m, 0:xt.get_shape()[1]] = xt
		xt_full.set_axisnam(["grp_"+str(i+1) for i in range(0, m)], 
							["obs_"+str(i+1) for i in range(0, xt_full.get_shape()[1])])
		design["xt"] = xt_full
		xt = design["xt"]
		design_new["xt"] = xt

		# model switch has to increase
		model_switch_full = Matrix(np.ones(m, int(np.amax(maxni.get_data()))) * np.nan)
		model_switch_full.get_data()[0:m, 0:model_switch.get_shape()[1]] = model_switch
		model_switch_full.set_axisnam(["grp_"+str(i+1) for i in range(0, m)],
									  ["obs_"+str(i+1) for i in range(0, model_switch_full.get_shape()[1])])
		design["model_switch"] = model_switch_full
		for i in range(0, design["model_switch"].get_shape()[0]):
			x_tmp = design["model_switch"].get_data().iloc[i, :]
			_, idx = np.unique(x_tmp[~np.isnan(x_tmp)], return_index=True) # remove duplicated and nan values but keep order
			if len([x_tmp[~np.isnan(x_tmp)][index] for index in sorted(idx)]) == 1:
				x_tmp[np.isnan(x_tmp)] = [x_tmp[~np.isnan(x_tmp)][index] for index in sorted(idx)]
			else:
				raise Exception("Unable to determine the model_switch values needed for group " + str(i+1)
								+ "\n Please supply them as input.")
			design["model_switch"].get_data()[i, :] = x_tmp
		model_switch = design["model_switch"]
		design_new["model_switch"] = model_switch


	# maxgroupsize
	if maxgroupsize.get_shape()[0] == 1 and m != 1:
		maxgroupsize.repeat([1, m], [m, 1], axisnam=[["grp_"+str(i+1) for i in range(0, m)], None])
	if test_mat_size([m, 1], maxgroupsize.get_shape(), "maxgroupsize") == 1:
		maxgroupsize.set_axisnam([["grp_"+str(i) for i in range(1, m+1)],
								["n_id"] * maxgroupsize.get_shape()[1]])
	
	# mingroupsize
	if mingroupsize.get_shape()[0] == 1 and m != 1:
		mingroupsize.repeat([1, m], [m, 1], axisnam=[["grp_"+str(i+1) for i in range(0, m)], None])
	if test_mat_size([m, 1], mingroupsize.get_shape(), "mingroupsize") == 1:
		mingroupsize.set_axisnam([["grp_"+str(i) for i in range(1, m+1)],
								["n_id"] * mingroupsize.get_shape()[1]])

	# make sure min is less than max
	ret = comp_max_min(maxgroupsize.get_data(), mingroupsize.get_data(), called_args)
	maxgroupsize = Matrix(ret["max_val"])
	mingroupsize = Matrix(ret["min_val"])

	# check given max and min
	if any(groupsize.get_data() < mingroupsize.get_data()):
		raise Exception("groupsize is less than mingroupsize")
	if any(groupsize.get_data() > maxgroupsize.get_data()):
		raise Exception("groupsize is greater than maxgroupsize")

	# maxtotgroupsize
	if maxtotgroupsize is None:
		maxtotgroupsize = Matrix(np.sum(groupsize.get_data()))
	
	# mintotgroupsize
	if mintotgroupsize is None:
		mintotgroupsize = Matrix(np.sum(mingroupsize.get_data()))

	# make sure min is less than max
	ret = comp_max_min(maxtotgroupsize.get_data(), mintotgroupsize.get_data(), called_args)
	maxtotgroupsize = Matrix(ret["max_val"])
	mintotgroupsize = Matrix(ret["min_val"])

	# check given max and min
	if any(np.sum("groupsize".get_data()) < mintotgroupsize.get_data()):
		raise Exception("sum of groupsizes is less than mintotgroupsize")
	if any(np.sum("groupsize".get_data()) > maxtotgroupsize.get_data()):
		raise Exception("sum of groupsizes is greater than maxtotgroupsize")
	
	# maxxt and minxt
	if type(maxxt) is int or type(maxxt) is float:
		maxxt = Matrix(np.ones(xt.get_shape()) * maxxt)
	if type(maxxt) is list:
		# length = max([len(i) for i in maxxt])
		# maxxt_ = []
		# for i in range(0, len(maxxt)):
		# 	maxxt[i] = maxxt[i].astype(np.float32)
		# 	maxxt[i] = np.pad(maxxt[i], (0, length-len(maxxt[i])), "constant", constant_values=np.nan)
		# 	maxxt_.append(maxxt[i].tolist())
		maxxt = Matrix(maxxt)
	if type(maxxt) is not Matrix:
		maxxt = Matrix(maxxt)
	if maxxt.get_shape()[0] == 1 and m != 1:
		maxxt.repeat([1, m], [m, maxxt.get_shape()[1]])
	if maxxt.get_shape()[0] != m:
		raise Exception("The number of rows in maxxt (" +
						str(maxxt.shape()[0]) +
						") is not the same as the number of groups m (" +
						str(m) +
						")")
	if maxxt.get_shape()[1] == int(np.max(ni.get_data())) and int(np.max(maxni.get_data())) > int(np.max(ni.get_data())) and xt.get_shape()[1] == int(np.max(maxni.get_data())):
		maxxt_full = xt.get_data()
		maxxt_full[:, 0:np.max(ni.get_data())] = maxxt.get_data()
		maxxt = Matrix(maxxt_full)
	if test_mat_size(xt.get_shape(), maxxt.get_shape(), "maxxt") == 1:
		maxxt.set_axisnam([["grp_"+str(i+1) for i in range(0, design["m"])], 
							["obs_"+str(i+1) for i in range(0, maxxt.shape[1])]])

	if type(minxt) is int or type(minxt) is float:
		minxt = Matrix(np.ones(xt.get_shape()) * minxt)
	if type(minxt) is list:
		# length = max([len(i) for i in minxt])
		# minxt_ = []
		# for i in range(0, len(minxt)):
		# 	minxt[i] = minxt[i].astype(np.float32)
		# 	minxt[i] = np.pad(minxt[i], (0, length-len(minxt[i])), "constant", constant_values=np.nan)
		# 	minxt_.append(minxt[i].tolist())
		minxt = Matrix(minxt)
	if type(minxt) is not Matrix:
		minxt = Matrix(minxt)
	if minxt.get_shape()[0] == 1 and m != 1:
		minxt.repeat([1, m], [m, minxt.get_shape()[1]])
	if minxt.get_shape()[0] != m:
		raise Exception("The number of rows in maxxt (" +
						str(minxt.shape()[0]) +
						") is not the same as the number of groups m (" +
						str(m) +
						")")
	if minxt.get_shape()[1] == int(np.max(ni.get_data())) and int(np.max(maxni.get_data())) > int(np.max(ni.get_data())) and xt.get_shape()[1] == int(np.max(maxni.get_data())):
		minxt_full = xt.get_data()
		minxt_full[:, 0:np.max(ni.get_data())] = minxt.get_data()
		minxt = Matrix(minxt_full)
	if test_mat_size(xt.get_shape(), minxt.get_shape(), "maxxt") == 1:
		minxt.set_axisnam([["grp_"+str(i+1) for i in range(0, design["m"])], 
							["obs_"+str(i+1) for i in range(0, minxt.shape[1])]])

	# make sure min is less than max
	ret = comp_max_min(maxxt, minxt, called_args)
	maxxt = Matrix(ret["max_val"])
	minxt = Matrix(ret["min_val"])

	# check for zeros
	if our_zero is not None:
		minxt = Matrix(minxt.get_data() + our_zero * (minxt.get_data() == 0))
		maxxt = Matrix(maxxt.get_data() + our_zero * (maxxt.get_data() == 0))
		design["xt"] = Matrix(xt.get_data() + our_zero * (xt.get_data() == 0))
		xt = design["xt"]
	
	# check given max and min
	if np.greater(minxt.get_data(), xt.get_data()).any():
		raise Exception("xt is less than minxt")
	if np.greater(xt.get_data(), maxxt.get_data()).any():
		raise Exception("xt is greater than maxxt")
	
	# need to decide on appripriate values of xt and minxt and maxxt if applicable
	if any(maxni.get_data() > ni.get_data()) and any(np.isnan(xt.get_data())):
		for grp in range(0, m):
			xt_grp = xt.get_data()[grp, :]
			maxxt_grp = maxxt.get_data()[grp, :]
			minxt_grp = minxt.get_data()[grp, :]
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
			design["xt"].get_data()[grp, :] = xt_grp
			xt = design["xt"]
			maxxt.get_data()[grp, :] = maxxt_grp
			minxt.get_data()[grp, :] = minxt_grp
		design_new["xt"] = xt

	# for a
	if maxa is not None:
		if type(maxa) is not Matrix:
			maxa = Matrix(maxa)
		if maxa.get_shape()[0] == 1 and m != 1:
			maxa.repeat([1, m], [m, maxa.get_size()])
		if maxa.get_shape()[0] != m:
			raise Exception("The number of rows in maxa (" +
							str(maxa.get_shape()[0]) +
							") is not the same as the number of groups m (" +
							str(m) + ")")
		colnam = maxa.get_axisnam()[1]
		if colnam is None:
			colnam = a.get_axisnam()[1]
		maxa.set_axisnam([["grp_" + str(i+1) for i in range(0, m)], colnam])
		design_space["maxa"] = maxa
	
	if mina is not None:
		if type(mina) is not Matrix:
			mina = Matrix(mina)
		if mina.get_shape()[0] == 1 and m != 1:
			mina.repeat([1, m], [m, mina.get_size()])
		if mina.get_shape()[0] != m:
			raise Exception("The number of rows in mina (" +
							str(mina.get_shape()[0]) +
							") is not the same as the number of groups m (" +
							str(m) + ")")
		colnam = mina.get_axisnam()[1]
		if colnam is None:
			colnam = a.get_axisnam()[1]
		mina.set_axisnam([["grp_" + str(i+1) for i in range(0, m)], colnam])
		design_space["mina"] = mina
	

	# make sure max is min smaller than max
	if mina is not None and maxa is not None:
		ret = comp_max_min(maxa, mina, called_args)
		maxa = Matrix(ret["max_val"])
		mina = Matrix(ret["min_val"])

	# check ni given max and min
	if mina is not None and maxa is not None and a is not None:
		if (np.greater(mina.get_data(), a.get_data())).any():
			raise Exception("a is less than mina")
		if (np.greater(a.get_data(), maxa.get_data())).any():
			raise Exception("a is greater than maxa")

	# for x
	if x_space is None and x is not None:
		x_space = Matrix(x.get_data()) # !!! not sure if it's right
	if x_space is not None:
		if x_space.get_shape()[0] == 1 and m != 1:
			x_space.repeat([1, m], [m, x_space.get_size()])
		if test_mat_size(x.get_shape(), x_space.get_shape(), "x_space") == 1:
			x_space.set_axisnam([["grp_"+str(i+1) for i in range(0, design["m"])],
									x.get_axisnam()[1]])
		design_space["x_space"] = x_space

		for i in range(0, x.get_shape()[0]):
			for j in range(0, x.get_shape()[1]):
				if type(x_space.get_one_data(index=[i, j])) is Matrix:
					tmp = x_space.get_one_data(index=[i, j])
				else:
					tmp = [x_space.get_one_data(index=[i, j])]
				if x.get_one_data(index=[i, j]) not in tmp:
					raise Exception("x value for group " + str(i+1) + " (column " + str(j+1) + ") is not in the design space")
	
	# for xt_space
	if xt_space is not None:
		if type(xt_space) is list:  # then it is a list with no dimensions, need to convert to a cell
			nspace = len(xt_space)
			nrow_xt = xt.get_shape()[0]
			ncol_xt = xt.get_shape()[1]
			if nspace == 1: # all time points in all groups have the same space
				xt_space = Matrix(xt_space).expand(xt.get_shape())
			elif nspace == ncol_xt: # we assume that all groups have the same space
				# len(xt_space) == xt.get_shape()[1] (column)
				xt_space_tmp = np.array(xt_space)
				xt_space = Matrix(xt_space).expand(xt.get_shape(), fill_value=0) # expand xt_space to xt size
				for i in range(0, nrow_xt):
					xt_space.get_data()[i, :] = xt_space_tmp # every new xt_space row is whole old xt_space values
			elif nspace == (ncol_xt * nrow_xt): # we assume that spaces are entered in row major form
				xt_space = Matrix(xt_space)
				xt_space.set_shape([int(xt_space.get_size()/nrow_xt), nrow_xt])
		else: # then it is a vector, assume the vector is the same for all xt's
			if type(xt_space) is not Matrix:
				xt_space = Matrix(xt_space)
			xt_space.expand(xt.get_shape())

		if xt_space.get_shape()[0] == 1 and m != 1:
			xt_space.repeat([1, m], [m, xt_space.get_size()])
		if xt_space.get_shape()[1] == 1 and xt.get_shape()[2] != 1:
			xt_space.repeat([xt.get_shape()[1], 1], [m, xt.get_shape()[1]])
		if test_mat_size(xt.get_shape(), xt_space.get_shape(0), "xt_space") == 1:
			xt_space.set_axisnam([["grp_"+str(i+1) for i in range(0, m)], xt.get_axisnam()[1]])
		
		for i in range(0, xt.get_shape()[0]):
			for j in range(0, xt.get_shape()[1]):
				if ~np.isnan(xt.get_data())[i, j]:
					if type(xt_space.get_data()[i, j]) is int or type(xt_space.get_data()[i, j]) is np.float64:
						tmp = [xt_space.get_data()[i, j]]
					else: # xt_space.get_data()[i, j] is list
						tmp = xt_space.get_data()[i, j]
					if xt.get_data()[i, j] not in tmp:
						raise Exception("xt value for group " + str(i+1) + " (column " + str(j+1) + ") is not in the design space")

	# for a_space
	if a_space is not None:
		if type(a_space) is list:
			if type(a_space[0]) is not list:
				a_space = Matrix(a_space * m, [m, len(a_space)])
			else:
				tmp_lst = []
				for i in range(0, len(a_space)):
					tmp_lst[i] = Matrix(a_space[i], [1, len(a_space[i])]) # !!! did not write "dimnames = list(None,names(x))"							
				mat = None
				for jj in range(0, len(tmp_lst)):
					tmp = tmp_lst[jj]
					if tmp.get_axisnam()[1] is not None and a.get_axisnam()[1] is not None:
						tmp = tmp.get_data()[:, [int(i) for i in a.get_axisnam()[1]]]
						mat = Matrix(tmp)
				a_space = mat

		if a_space.get_shape()[0] == 1 and m != 1:
			a_space.repeat([1, m], [m, a_space.get_size()])
		if a_space.get_shape()[1] == 1 and a.get_shape()[1] != 1:
			a_space.repeat([a.get_shape()[1], 1], [m, a.get_shape(0[1])])
		if test_mat_size(a.get_shape(), a_space.get_shape(), "a_space") == 1:
			a_space.set_axisnam([["grp_" + str(i+1) for i in range(0, m)], a.get_axisnam()[1]])
		for i in range(0, a.get_shape()[0]):
			for j in range(0, a.get_shape()[1]):
				if a.get_data()[i, j] is not None:
					if a.get_data()[i, j] not in [a_space.get_data()[i, j]]
						raise Exception("a value for group " + str(i+1) + " (column " + str(j+1) + ") is not in the design space")

	# for grouped_xt
	if grouped_xt is None:
		grouped_xt = Matrix(design["xt"].get_data() * np.nan)
		val = 1
		for i in range(0, design["xt"].get_shape()[0]):
			if use_grouped_xt:
				val = 1
			for j in range(0, design["xt"].get_shape()[1]):
				if design["xt"].get_data()[i, j] is not None:
					grouped_xt.set_one_data(val, index=[i, j])
					val += 1

	if grouped_xt.size == 1:
		grouped_xt = Matrix(np.array(ones(size(design["xt"])[0], size(design["xt"])[1])) * grouped_xt)
		use_grouped_xt = True
	if type(grouped_xt) is list:
		length = max([len(i) for i in grouped_xt])
		grouped_xt_ = []
		for i in range(0, len(grouped_xt)):
			grouped_xt[i] = grouped_xt[i].astype(np.float32)
			grouped_xt[i] = np.pad(grouped_xt[i], (0, length-len(grouped_xt[i])), "constant", constant_values=np.nan)
			grouped_xt_.append(grouped_xt[i].tolist())
		grouped_xt = Matrix(grouped_xt_)
	if size(grouped_xt)[0] == 1 and design["m"] != 1:
		grouped_xt = Matrix(np.tile(grouped_xt.flatten(), design["m"]).reshape(design["m"], grouped_xt.size))
		use_grouped_xt = True
	if type(grouped_xt) is not Matrix:
		grouped_xt = Matrix(grouped_xt)
	if size(grouped_xt)[1] == np.max(design["ni"].get_data()) and np.max(maxni.get_data()) > np.max(design["ni"].get_data()) and size(design["xt"])[1] == np.max(maxni.get_data()):
		grouped_xt_full = design["xt"].get_data() * np.nan
		grouped_xt_full[:, 0:int(np.max(design["ni"]))] = np.array(grouped_xt.astype(np.float32))
		grouped_xt = Matrix(grouped_xt_full)
	if test_mat_size(np.array(size(design["xt"])), grouped_xt.get_data(), "grouped_xt") == 1:
		grouped_xt = Matrix(grouped_xt,
								  rownam=["grp_" + str(i+1) for i in range(0, design["m"])],
								  colnam=["obs_" + str(i+1) for i in range(0, grouped_xt.shape[1])])

	# get values in the NA region if possible
	if any(maxni.get_data() > design["ni"].get_data()) and any(np.isnan(grouped_xt.get_data())):
		for grp in range(0, design["m"]):
			grouped_xt_grp = grouped_xt.get_data()[grp, :]
			if any(np.isnan(grouped_xt_grp)):
				_, idx = np.unique(grouped_xt_grp[~np.isnan(grouped_xt_grp)], return_index=True) # remove duplicated and nan values but keep order
				vals = [grouped_xt_grp[~np.isnan(grouped_xt_grp)][index] for index in sorted(idx)]
				if len(vals) == 1:
					grouped_xt_grp[np.isnan(grouped_xt_grp)] = vals
				else:
					raise Exception("Unable to determine the grouped_xt values needed for group " +
									str(grp+1) +
									"\n if ni increases with optimization \nPlease supply them as input.")
			grouped_xt.set_multiple_data(grouped_xt_grp, row=grp)

	_, idx = np.unique(grouped_xt.get_data()[~np.isnan(design["xt"].get_data())], return_index=True) # remove duplicated and nan values but keep order
	for i in [grouped_xt.get_data()[~np.isnan(design["xt"].get_data())][index] for index in sorted(idx)]:
		_, idx = np.unique(design["xt"].get_data()[np.logical_and(grouped_xt.get_data() == i, ~np.isnan(design["xt"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
		if len([design["xt"].get_data()[np.logical_and(grouped_xt.get_data() == i, ~np.isnan(design["xt"].get_data()))][index] for index in sorted(idx)]) != 1:
			raise Exception("xt values grouped with value %g from grouped_xt do not have the same initial values.\n'" % i)
		_, idx = np.unique(maxxt.get_data()[np.logical_and(grouped_xt.get_data() == i, ~np.isnan(design["xt"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
		if len([maxxt.get_data()[np.logical_and(grouped_xt.get_data() == i, ~np.isnan(design["xt"].get_data()))][index] for index in sorted(idx)]) != 1:
			raise Exception("xt values grouped with value %g from grouped_xt do not have the same maximum allowed values (maxxt).\n" % i)
		_, idx = np.unique(minxt.get_data()[np.logical_and(grouped_xt.get_data() == i, ~np.isnan(design["xt"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
		if len([minxt.get_data()[np.logical_and(grouped_xt.get_data() == i, ~np.isnan(design["xt"].get_data()))][index] for index in sorted(idx)]) != 1:
			raise Exception("xt values grouped with value %g from grouped_xt do not have the same maximum allowed values (minxt).\n" % i)
		if xt_space is None:
			grouped_cells_xt = None
		else:
			grouped_cells_xt = xt_space.get_data()[np.logical_and(grouped_xt.get_data() == i, ~np.isnan(design["xt"].get_data()))]
			for j in range(0, grouped_cells_xt.get_shape()[0]):
				for k in range(j, grouped_cells_xt.get_shape()[0]):
					if ((np.array(size(grouped_cells_xt.get_data()[j])) != np.array(size(grouped_cells_xt.get_data()[k]))).any() or
						(grouped_cells_xt.get_data()[j] != grouped_cells_xt.get_data()[k]).any()):
						raise Exception("xt values grouped with value % g from grouped_xt do not have the same allowed discrete values (xt_space).\n" % i)

	_, idx = np.unique(grouped_xt.get_data()[~np.isnan(design["xt"].get_data())], return_index=True) # remove duplicated and nan values but keep order
	for i in range(0, int(np.max([grouped_xt.get_data()[~np.isnan(design["xt"].get_data())][index] for index in sorted(idx)]))):
		_, idx = np.unique(design["xt"].get_data()[np.logical_and(grouped_xt.get_data() == i+1, ~np.isnan(design["xt"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
		if len([design["xt"].get_data()[np.logical_and(grouped_xt.get_data() == i+1, ~np.isnan(design["xt"].get_data()))][index] for index in sorted(idx)]) == 0:
			raise Exception("grouped_xt must be sequential and cannot have missing values.\nNo xt values were grouped with value %g in grouped_xt.\n" % i)

	# for grouped_a
	if "a" in design.keys():
		if design["a"] is not None:
			if grouped_a is None:
				grouped_a = Matrix(design["a"].get_data() * np.nan)
				val = 1
				for i in range(0, size(design["a"])[0]):
					if use_grouped_a:
						val = 1
					for j in range(0, size(design["a"])[1]):
						if ~np.isnan(design["a"].get_one_data(index=[i, j])):
							grouped_a.set_one_data(val, index=[i, j])
							val += 1

			if size(grouped_a) == [1, 1]:
				grouped_a = Matrix(np.ones(size(design["a"])) * grouped_a.get_one_data(index=[0, 0]))

			if type(grouped_a) == list:
				length = max([len(i) for i in grouped_a])
				grouped_a_ = []
				for i in range(0, len(grouped_a)):
					grouped_a[i] = grouped_a[i].astype(np.float32)
					grouped_a[i] = np.pad(grouped_a[i], (0,length-len(grouped_a[i])),"constant",constant_values=np.nan)
					grouped_a_.append(grouped_a[i].tolist())
				grouped_a = Matrix(np.array(grouped_a_))
				use_grouped_a = True

			if type(grouped_a) is not Matrix:
				grouped_a = Matrix(grouped_a)

			if test_mat_size(np.array(size(design["a"])), grouped_a, "grouped_a") == 1:
				grouped_a = Matrix(grouped_a,
										rownam=["grp_" + str(i+1) for i in range(0, design["m"])])
				if grouped_a.get_colnam == []:
					grouped_a.set_colnam(design["a"].get_colnam())

			_, idx = np.unique(grouped_a.get_data()[~np.isnan(design["a"].get_data())], return_index=True) # remove duplicated and nan values but keep order
			for i in [grouped_a.get_data()[~np.isnan(design["a"].get_data())][index] for index in sorted(idx)]:
				_, idx = np.unique(design["a"].get_data().reshape(design["a"].get_shape())[np.logical_and(grouped_a.get_data() == i, ~np.isnan(design["a"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
				if len([design["a"].get_data().reshape(design["a"].get_shape())[np.logical_and(grouped_a.get_data() == i, ~np.isnan(design["a"].get_data()))][index] for index in sorted(idx)]) != 1:
					raise Exception("a values grouped with value %g from grouped_a do not have the same initial values.\n'" % i)
				_, idx = np.unique(maxa.get_data().reshape(maxa.get_shape())[np.logical_and(grouped_a.get_data() == i, ~np.isnan(design["a"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
				if len([maxa.get_data().reshape(maxa.get_shape())[np.logical_and(grouped_a.get_data() == i, ~np.isnan(design["a"].get_data()))][index] for index in sorted(idx)]) != 1:
					raise Exception("a values grouped with value %g from grouped_a do not have the same maximum allowed values (maxa).\n" % i)
				_, idx = np.unique(mina.get_data().reshape(mina.get_shape())[np.logical_and(grouped_a.get_data() == i, ~np.isnan(design["a"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
				if len([mina.get_data().reshape(mina.get_shape())[np.logical_and(grouped_a.get_data() == i, ~np.isnan(design["a"].get_data()))][index] for index in sorted(idx)]) != 1:
					raise Exception("a values grouped with value %g from grouped_a do not have the same maximum allowed values (mina).\n" % i)
				if a_space is None:
					grouped_cells_a = None
				else:
					grouped_cells_a = Matrix(a_space.get_data()[np.logical_and(grouped_a.get_data() == i, np.isnan(design["a"].get_data()))])
					for j in range(0, grouped_cells_a.get_size()):
						for k in range(0, grouped_cells_a.get_size()):
							if (any(np.array(size(grouped_cells_a[j])) != np.array(size(grouped_cells_a[k]))) or
								(grouped_cells_a.get_data()[j] != grouped_cells_a.get_data()[k])): # 没写any
								raise Exception("a values grouped with value %g from grouped_a do not have the same allowed discrete values (a_space).\n" % i)

			_, idx = np.unique(grouped_a.get_data()[~np.isnan(design["a"].get_data())], return_index=True) # remove duplicated and nan values but keep order
			for i in range(0, int(np.max([grouped_a.get_data()[~np.isnan(design["a"].get_data())][index] for index in sorted(idx)]))):
				_, idx = np.unique(design["a"].get_data().reshape(design["a"].get_shape())[np.logical_and(grouped_a.get_data() == i+1, ~np.isnan(design["a"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
				if len([design["a"].get_data().reshape(design["a"].get_shape())[np.logical_and(grouped_a.get_data() == i+1, ~np.isnan(design["a"].get_data()))][index] for index in sorted(idx)]) == 0:
					raise Exception("grouped_a must be sequential and cannot have missing values.\nNo a values were grouped with value %g in grouped_a.\n" % i)

			design_space["grouped_a"] = grouped_a
			design_space["use_grouped_a"] = use_grouped_a

	# for grouped_x
	if "x" in design.keys():
		if design["x"] is not None:
			if grouped_x is None:
				grouped_x = Matrix(np.array(design["x"]) * np.nan)
				val = 1
				for i in range(0, size(design["x"])[0]):
					if use_grouped_x:
						val = 1
					for j in range(0, size(design["x"])[1]):
						if ~np.isnan(design["x"].get_data()[i, j]):
							grouped_x[i, j] = val
							val += 1

			if len(grouped_x) == 1:
				grouped_x = Matrix(np.array(ones(size(design["x"].get_data()[0]), size(design["x"][1]))) * grouped_x)

			if type(grouped_x) == list:
				length = max([len(i) for i in grouped_x])
				grouped_x_ = []
				for i in range(0, len(grouped_x)):
					grouped_x[i] = grouped_x[i].astype(np.float32)
					grouped_x[i] = np.pad(grouped_x[i], (0,length-len(grouped_x[i])),"constant",constant_values=np.nan)
					grouped_x_.append(grouped_x[i].tolist())
				grouped_x = Matrix(np.array(grouped_x_))
				use_grouped_x = True

			if type(grouped_x) is not Matrix:
				grouped_x = Matrix(grouped_x)

			if test_mat_size(np.array(size(design["x"])), grouped_x.get_data(), "grouped_x") == 1:
				grouped_x.set_rownam(["grp_" + str(i+1) for i in range(0, design["m"])])
				if grouped_x.get_colnam() == []:
					grouped_x.set_colnam(design["x"].get_colnam())

			_, idx = np.unique(grouped_x.get_data()[~np.isnan(design["x"].get_data())], return_index=True) # remove duplicated and nan values but keep order
			for i in [grouped_x.get_data()[~np.isnan(design["x"].get_data())][index] for index in sorted(idx)]:
				_, idx = np.unique(design["x"].get_data()[np.logical_and(grouped_x.get_data() == i, ~np.isnan(design["x"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
				if len([design["x"][np.logical_and(grouped_x.get_data() == i, ~np.isnan(design["x"].get_data()))][index] for index in sorted(idx)]) != 1:
					raise Exception("x values grouped with value %g from grouped_x do not have the same initial values.\n'" % i)
				grouped_cells = x_space.get_data()[np.logical_and(grouped_x.get_data() == i, np.isnan(design["x"].get_data()))]
				for j in range(0, grouped_cells.size):
					for k in range(0, grouped_cells.size):
						if (any(np.array(size(grouped_cells[j])) != np.array(size(grouped_cells[k]))) or
							grouped_cells[j] != grouped_cells[k]):
							raise Exception("x values grouped with value %g from grouped_x do not have the same allowed discrete values (grouped_x).\n" % i)

			_, idx = np.unique(grouped_x.get_data()[~np.isnan(design["x"].get_data())], return_index=True) # remove duplicated and nan values but keep order
			for i in range(0, int(np.max([grouped_x.get_data()[~np.isnan(design["x"].get_data())][index] for index in sorted(idx)]))):
				_, idx = np.unique(design["x"].get_data()[np.logical_and(grouped_x.get_data() == i+1, ~np.isnan(design["x"].get_data()))], return_index=True) # remove duplicated and nan values but keep order
				if len([design["x"].get_data()[np.logical_and(grouped_x.get_data() == i+1, ~np.isnan(design["x"].get_data()))][index] for index in sorted(idx)]) == 0:
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
		for i in range(0, a_space.get_shape()[0]):
			for j in range(0, a_space.get_shape()[1]):
				maxa[i, j] = np.max(a_space.get_data()[i, j][0])
	if mina_imputed and a_space is not None:
		for i in range(0, a_space.get_shape()[0]):
			for j in range(0, a_space.get_shape()[1]):
				mina[i, j] = np.min(a_space.get_data()[i, j][0])
	design_space["maxa"] = maxa
	design_space["mina"] = mina

	if maxxt_imputed and xt_space is not None:
		for i in range(0, xt_space.get_shape()[0]):
			for j in range(0, xt_space.get_shape()[1]):
				maxxt[i, j] = np.max(xt_space.get_data()[i, j][0])
	if minxt_imputed and xt_space is not None:
		for i in range(0, xt_space.get_shape()[0]):
			for j in range(0, xt_space.get_shape()[1]):
				minxt[i, j] = np.min(xt_space.get_data()[i, j][0])
	design_space["maxxt"] = maxxt
	design_space["minxt"] = minxt

	return {"design": design_new, "design_space": design_space}

def comp_max_min(max_val, min_val, called_args) -> np.ndarray:
	args = list(locals().values()) # []: argument values
	if type(min_val) is Matrix and type(max_val) is Matrix:
		if np.greater(min_val.get_data(), max_val.get_data()).any(): # if any of min_val is greater than any of max_val
			# check if args[0] and args[1] in called_args key
			min_val_sup = str(args[1]) in list(called_args.keys())
			max_val_sup = str(args[0]) in list(called_args.keys())
			if min_val_sup and max_val_sup:
				raise Exception("Some value of " + str(args[0]) + " is smaller than " + args[1])
			if min_val_sup and ~max_val_sup:
				max_val = np.maximum(max_val.get_data(), min_val.get_data())
			if ~min_val_sup and max_val_sup:
				min_val = np.minimum(max_val.get_data(), min_val.get_data())
	elif type(min_val) is np.ndarray and type(max_val) is np.ndarray:
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