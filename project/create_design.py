"""
Author: Caiya Zhang, Yuchen Zheng
"""


#import project.all_modules as am
from enum import IntFlag
import re
import numpy as np
import pandas as pd
import itertools
from matpy.matrix import matrix
from project.poped_choose import poped_choose
from project.test_mat_size import test_mat_size
from project.size import size


def create_design(
				  xt, # Matrix defining the initial sampling schedule row major, can also be a list of vectors
				  groupsize, # Vector defining the size of the different groups (num individuals in each group)
				  m=None, # Number of groups, computed from xt if not defined
				  x=None, # Matrix defining the initial discrete values
				  a=None,
				  ni=None, # Vector defining the number of samples for each group, computed as all elements of xt by default
				  model_switch=None # Vector defining which response a certain sampling time belongs to, defaults to one for all elements of xt
):

	design = {}

	### for xt, m ###
	if type(xt) is list: 
		length = max([len(i) for i in xt])
		xt_ = []
		for i in range(0, len(xt)):
			xt[i] = xt[i].astype(np.float32)
			xt[i] = np.pad(xt[i], (0, length-len(xt[i])), "constant", constant_values=np.nan)
			# xt = np.array([np.pad(i, (0, length-len(i)), 'constant', constant_values=np.nan) for i in xt]) # convert a list of vectors to an array
			xt_.append(xt[i].tolist())
		xt = matrix(np.array(xt_))

	if m is None: 
		m = xt.get_shape()[0] # get xt row (same as "m = size(xt, 1)")

	if size(xt)[0] == 1 and m != 1:
		xt = matrix(np.tile(xt.get_all_data().flatten(), m), shape=(m, xt.size)) # flatten xt, repeat by m times, and reshape to (col: xt's element number, row: m)

	if type(xt) is not matrix:
		xt = matrix(xt)
	
	if (size(xt)[0] != m):
		raise Exception("The number of rows in xt (" + str(size(xt)[0]) + ") is not the same as the number of groups m (" + str(m) + ")")
	
	xt.set_rownam(["grp_"+str(i) for i in range(1, m+1)]) 
	xt.set_colnam(["obs_"+str(i) for i in range(1, xt.shape[1]+1)]) # same as "size(xt)[1]+1"
	
	design["xt"] = xt
# 没写！！！names(m) = "n_grp"
	design["m"] = m


	### for ni ###
	if ni is None:
		ni = np.count_nonzero(1-np.isnan(xt.get_all_data()), axis=1).reshape(xt.get_shape()[0], 1)
	
	if type(ni) is not matrix:
		ni = matrix(np.array(ni), shape=[len(ni), 1])
	if test_mat_size(np.array([m, 1]), ni.get_all_data(), "ni") == 1:
		ni.set_axisnam(list(itertools.repeat("n_obs", ni.shape[1])),
					   ["grp_"+str(i) for i in range(1, m+1)]) 
		design["ni"] = ni


	### for model_switch ###
	if type(model_switch) is list:
		length = max([len(i) for i in model_switch])
		model_switch = matrix(np.array([np.pad(i, (0, length-len(i)), 'constant', constant_values=np.nan) for i in model_switch])) # convert a list of vectors to an array
	if model_switch is None:
		model_switch = matrix(xt.get_all_data() * 0 + 1)
	if size(model_switch)[0] == 1 and m != 1:
		model_switch  = matrix(np.tile(model_switch.get_all_data().flatten(), m), shape=(m, model_switch.size)) # flatten xt, repeat by m times, and reshape to (col: xt's element number, row: m)

	if type(model_switch) is not matrix:
		model_switch = matrix(model_switch)

	if test_mat_size(np.array(size(xt)), model_switch.get_all_data(), "model_switch") == 1:
		model_switch.set_axisnam(["obs_"+str(i) for i in range(1, model_switch.shape[1]+1)],
								 ["grp_"+str(i) for i in range(1, m+1)]) 
		design["model_switch"] = model_switch


	### for a ###
	if a is not None:
		if type(a) is list or type(a) is not matrix:
			a = matrix(a)
		# elif len(size(a)) == 1:
		# 	a = np.array([a])
		colnam = None
		if size(a)[0] == 1 and m != 1:
			a_ = []
			if type(a) is int or size(a) == [1, 1]:
				a = matrix(np.tile([a], m), shape=(m, 1),
						   colnam=colnam, rownam=["grp_"+str(i) for i in range(1, m+1)])
			else:
				for i in range(0, size(a)[1]):
					for j in range(0, size(a)[0]):
						a_.append(a.get_all_data()[j][i])
				a = matrix(np.tile(a_, m), shape=(m, a.size),
						   colnam=colnam, rownam=["grp_"+str(i) for i in range(1, m+1)])

		a = matrix(a)
		if size(a)[0] != m:
			raise Exception("The number of rows in a (" + str(size(a)[0]) + ") is not the same as the number of groups m (" + str(m) + ")")
		a.set_rownam(["grp_"+str(i) for i in range(1, m+1)])
		if a.get_colnam() is not None:
			count = 0
			for i in range(0, a.get_shape()[1]):
				if re.search("^X[0-9]+$", str(a.get_colnam()[i])) is not None:
					count += 1
			if count == size(a)[1]:
				a.set_colnam([None] * a.shape[1])
		design["a"] = a


	### for x ###
	if x is not None:
		if type(x) == list:
			x = matrix(x)
		colnam = x.get_colnam()
		if size(x)[0] == 1 and m != 1:
			x_ = []
			for i in range(0, x.get_shape()[1]):
				for j in range(0, x.get_shape()[0]):
					x_.append(x.get_all_data()[i][j])
			x = matrix(np.tile(x_, m), shape=(m, x.size),
					   colnam=colnam, rownam=["grp_"+str(i) for i in range(1, m+1)])

		if type(x) is not matrix:
			x = matrix(x)

		if size(x)[0] != m:
			raise Exception("The number of rows in x (" + str(size(x)[0]) + "is not the same as the number of groups m (" + str(m) + ")")
		x.set_rownam(["grp_"+str(i) for i in range(1, m+1)])
		design["x"] = x


	### for groupsize ###
	if max(size(groupsize)) == 1 and m != 1:
		groupsize = [groupsize] * m
		groupsize = matrix(np.array(groupsize), shape=[m ,1],
						   rownam=["grp_"+str(i) for i in range(1, m+1)])
	
	if type(groupsize) is matrix:
		if len(groupsize.get_shape()) != 2:
			groupsize = groupsize.get_all_data()[:, np.newaxis]
	if type(groupsize) is np.ndarray:
		if len(groupsize.get_shape()) != 2:
			groupsize = matrix(groupsize.get_all_data()[:, np.newaxis])
	elif type(groupsize) is list or type(groupsize) is int:
		groupsize = matrix(np.array([groupsize])[:, np.newaxis])
		
	if test_mat_size(np.array([m, 1]), groupsize.get_all_data(), "groupsize") == 1:
		groupsize = matrix(groupsize,
						   rownam=["grp_"+str(i) for i in range(1, m+1)],
						   colnam=["n_id"] * groupsize.shape[1])
		design["groupsize"] = groupsize
	

	return design