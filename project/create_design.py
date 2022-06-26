"""
Author: Caiya Zhang, Yuchen Zheng
"""


import path
from enum import IntFlag
import re
import numpy as np
import pandas as pd
import itertools
from matpy.matrix import Matrix, select
from matpy.num import Num
from project.test_mat_size import test_mat_size


def create_design(
				  xt: Matrix, # Matrix defining the initial sampling schedule row major
				  groupsize, # Vector defining the size of the different groups (num individuals in each group)
				  m=None, # Number of groups, computed from xt if not defined
				  x=None, # Matrix defining the initial discrete values
				  a=None,
				  ni=None, # Vector defining the number of samples for each group, computed as all elements of xt by default
				  model_switch=None # Vector defining which response a certain sampling time belongs to, defaults to one for all elements of xt
):

	design = {}

	### for xt, m ###
	if m is None: 
		m = xt.get_shape()[0]
	if xt.get_shape()[0] == 1 and m != 1:
		xt.repeat([1, m], [m, xt.get_size()], True, False)
	if (xt.get_shape()[0] != m): # check
		raise Exception("The number of rows in xt (" + str(xt.get_shape()[0]) + ") is not the same as the number of groups m (" + str(m) + ")")
	xt.set_axisnam([["grp_"+str(i) for i in range(1, m+1)], ["obs_"+str(i) for i in range(1, xt.get_shape()[1]+1)]]) 
	design["xt"] = xt
	m = Num(m, "n_grp")
	design["m"] = m


	### for ni ###
	if ni is None:
		tmp: np.ndarray = 1 - np.isnan(xt.get_data())
		if len(tmp.shape) == 1:
			tmp = tmp[np.newaxis, :]
		ni = Matrix(np.count_nonzero(tmp, axis=1), [xt.get_shape()[0], 1])
	if test_mat_size([m.get_value(), 1], ni.get_shape(), "ni") == 1:
		ni.set_axisnam([["grp_"+str(i) for i in range(1, m.get_value()+1)], 
						list(itertools.repeat("n_obs", ni.get_shape()[1]))]) 
		design["ni"] = ni


	### for model_switch ###
	if model_switch is None:
		model_switch = Matrix(xt.get_data() * 0 + 1)
	if model_switch.get_shape()[0] == 1 and m.get_value() != 1:
		model_switch.repeat([1, m.get_value()], [m.get_value(), model_switch.get_size()], True, False)
	if test_mat_size(xt.get_shape(), model_switch.get_shape(), "model_switch") == 1:
		model_switch.set_axisnam([["grp_"+str(i) for i in range(1, m.get_value()+1)], 
								["obs_"+str(i) for i in range(1, model_switch.get_shape()[1]+1)]]) 
		design["model_switch"] = model_switch


	### for a ###
	if a is not None:
		if type(a) is not Matrix:
			a = Matrix(a)
	if type(a) is Matrix:
		if a.get_shape()[0] == 1 and m.get_value() != 1:
			a.repeat([1, m.get_value()], [m.get_value(), a.get_size()], datanam=True, axisnam=True)
		if a.get_shape()[0] != m.get_value():
			raise Exception("The number of rows in a (" + str(a.get_shape()[0]) + ") is not the same as the number of groups m (" + str(m.get_value()) + ")")
		colnam = list(range(0, a.get_shape()[1]))
		if a.get_axisnam() is not None:
			if a.get_axisnam()[1] is not None:
				colnam = a.get_axisnam()[1]
		a.set_axisnam([["grp_"+str(i) for i in range(1, m.get_value()+1)], colnam])
		if a.get_axisnam()[1] is not None:
			count = 0
			for i in range(0, a.get_shape()[1]):
				if re.search("^X[0-9]+$", str(a.get_axisnam()[1][i])) is not None:
					count += 1
			if count == a.get_shape()[1]:
				a.set_axisnam([a.get_axisnam()[0], [None] * a.get_shape()[1]])
		design["a"] = a


	### for x ###
	if type(x) is Matrix:
		# use datanam to create axisnam (deleted)
		if x.get_shape()[0] == 1 and m.get_value() != 1:
			x.repeat([1, m.get_value()], [m.get_value(), x.get_size()], datanam=True)
		if x.get_shape()[0] != m.get_value():
			raise Exception("The number of rows in x (" + str(x.get_shape()[0]) + "is not the same as the number of groups m (" + str(m.get_value()) + ")")
		x.set_axisnam(["grp_"+str(i) for i in range(1, m.get_value()+1)], x.get_axisnam()[1])
		design["x"] = x


	### for groupsize ###
	if type(groupsize) is not Matrix:
		groupsize = Matrix(groupsize)
	if max(groupsize.get_shape()) == 1 and m.get_value() != 1:
		groupsize.repeat([m.get_value(), 1], [m.get_value(), 1], datanam=True)
	else:
		groupsize.set_shape([m.get_value(), 1])	
	if test_mat_size([m.get_value(), 1], groupsize.get_shape(), "groupsize") == 1:
		groupsize.set_axisnam([["grp_"+str(i) for i in range(1, m.get_value()+1)],
								["n_id"] * groupsize.shape[1]])
		design["groupsize"] = groupsize
	
	return design