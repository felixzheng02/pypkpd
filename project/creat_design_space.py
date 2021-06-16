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
import inspect
from project.size import size


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

	args = inspect.getargspec(create_design_space).args

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
	
	design_space = []
	design_new = design

	for key, value in design.items():
		# rules:
		# 1. set default if not already defined
		# 2. read in value and translate to correct format
		# 3. check that size of object is correct
		# 4. add row and column names
		# 4. check that max is greater than min
		# 5. check that design value is wihin range of design_space
		if size(maxni, 1) == 1 and m != 1:
			...


def comp_max_min(max_val, min_val, called_args):
	args = inspect.getargspec(comp_max_min).args
	if any(np.greater(min_val, max_val)):
		 min_val_sup = args[3] in called_args
		 ...