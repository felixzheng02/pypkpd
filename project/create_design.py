"""
Author: Caiya Zhang, Yuchen Zheng
"""


#import project.all_modules as am
from enum import IntFlag
import re
import numpy as np
import pandas as pd
import itertools
from matpy.matrix import Matrix
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
		m = xt.get_shape()[0] # get xt row (same as "m = size(xt, 1)")

	if xt.get_shape()[0] == 1 and m != 1:
		xt.repeat([1, m], [m, xt.get_size()], True, False) # flatten xt, repeat by m times, and reshape to (col: xt's element number, row: m)
	
	if (xt.get_shape()[0] != m):
		raise Exception("The number of rows in xt (" + str(xt.get_shape()[0]) + ") is not the same as the number of groups m (" + str(m) + ")")
	
	xt.set_axisnam([["grp_"+str(i) for i in range(1, m+1)], ["obs_"+str(i) for i in range(1, xt.get_shape()[1]+1)]]) 
	
	design["xt"] = xt
	m = Num(m, "n_grp")
	design["m"] = m


	### for ni ###
	if ni is None:
		tmp = 1-np.isnan(xt.get_all_data())
		if len(tmp.shape) == 1:
			tmp = tmp[np.newaxis, :]
		ni = Matrix(np.count_nonzero(tmp, axis=1).reshape(xt.get_shape()[0], 1))
	
	
	if test_mat_size(np.array([m, 1]), ni.get_all_data(), "ni") == 1:
		ni.set_axisnam([["grp_"+str(i) for i in range(1, m+1)], 
						list(itertools.repeat("n_obs", ni.get_shape()[1]))]) 
		design["ni"] = ni


	### for model_switch ###
	# if type(model_switch) is list:
	# 	length = max([len(i) for i in model_switch])
	# 	model_switch = matrix(np.array([np.pad(i, (0, length-len(i)), 'constant', constant_values=np.nan) for i in model_switch])) # convert a list of vectors to an array
	if model_switch is None:
		model_switch = Matrix(xt.get_all_data() * 0 + 1)
	if model_switch.get_shape()[0] == 1 and m != 1:
		model_switch.repeat([1, m], [m, model_switch.get_size()], True, False)  # flatten xt, repeat by m times, and reshape to (col: xt's element number, row: m)

	if test_mat_size(np.array(xt.get_shape()), model_switch.get_all_data(), "model_switch") == 1:
		model_switch.set_axisnam([["grp_"+str(i) for i in range(1, m+1)], 
								["obs_"+str(i) for i in range(1, model_switch.get_shape()[1]+1)]]) 
		design["model_switch"] = model_switch


	### for a ###
	if type(a) is Matrix: # a is a Matrix
		colnam = a.get_datanam()
		if colnam is None:
			colnam = a.get_axisnam()[1]
		if a.get_shape()[0] == 1 and m != 1:
			a.repeat([1, m], [m, a.get_size()], datanam=True)
			a.set_axisnam(([["grp_"+str(i) for i in range(1, m+1)],
							colnam]))

		if a.get_shape()[0] != m:
			raise Exception("The number of rows in a (" + str(a.get_shape()[0]) + ") is not the same as the number of groups m (" + str(m) + ")")
		a.set_axisnam([["grp_"+str(i) for i in range(1, m+1)], a.get_axisnam()[1]])
		if a.get_axisnam()[1] is not None:
			count = 0
			for i in range(0, a.get_shape()[1]):
				if re.search("^X[0-9]+$", str(a.get_axisnam()[1][i])) is not None:
					count += 1
			if count == a.get_shape()[1]:
				a.set_axisnam([a.get_axisnam()[0], [None] * a.shape[1]])
		design["a"] = a


	### for x ###
	if type(x) is Matrix:
		colnam = x.get_datanam()
		if colnam is None:
			colnam = x.get_axisnam()[1]
		if size(x)[0] == 1 and m != 1:
			x.repeat([1, m], [m, x.get_size()], datanam=True)
			x.set_axisnam([["grp_"+str(i) for i in range(1, m+1)], colnam])

		if x.get_shape()[0] != m:
			raise Exception("The number of rows in x (" + str(x.get_shape()[0]) + "is not the same as the number of groups m (" + str(m) + ")")
		x.set_axisnam(["grp_"+str(i) for i in range(1, m+1)], x.get_axisnam()[1])
		design["x"] = x


	### for groupsize ###
	if type(groupsize) is Matrix:
		if max(groupsize.get_shape()) == 1 and m != 1:
			groupsize.repeat([m, 1], [m, 1], datanam=True)
			groupsize.set_axisnam([["grp_"+str(i) for i in range(1, m+1)], 
									None])
		
	if test_mat_size(np.array([m, 1]), groupsize.get_all_data(), "groupsize") == 1:
		groupsize.set_axisnam([["grp_"+str(i) for i in range(1, m+1)],
								["n_id"] * groupsize.shape[1]])
		design["groupsize"] = groupsize
	
	return design