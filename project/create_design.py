"""


Author: Caiya Zhang, Yuchen Zheng
"""


#import project.all_modules as am
import re
import numpy as np
import pandas as pd
import itertools
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
		xt = np.array(xt_)

	if m is None: 
		m = xt.shape[0] # get xt row (same as "m = size(xt, 1)")

	if size(xt)[0] == 1 and m != 1:
		xt = np.tile(xt.flatten(), m).reshape(m, xt.size) # flatten xt, repeat by m times, and reshape to (col: xt's element number, row: m)

	if type(xt) is not np.ndarray and type(xt) is not pd.DataFrame:
		xt = np.array(xt)
	
	if (size(xt)[0] != m):
		raise Exception("The number of rows in xt (" + str(size(xt)[0]) + ") is not the same as the number of groups m (" + str(m) + ")")
	
	xt = pd.DataFrame(xt, 
					   index=["grp_"+str(i) for i in range(1, m+1)], 
					   columns=["obs_"+str(i) for i in range(1, xt.shape[1]+1)]) # same as "size(xt)[1]+1"
	
	design["xt"] = xt
# 没写！！！names(m) <- "n_grp"
	design["m"] = m


	### for ni ###
	if ni is None:
		ni = np.count_nonzero(1-np.isnan(xt), axis=1).reshape(xt.shape[0], 1)
	
	if type(ni) != np.ndarray:
		ni = np.array(ni).reshape([len(ni), 1])
	if test_mat_size(np.array([m, 1]), ni, "ni") == 1:
		ni = pd.DataFrame(ni, 
					   index=["grp_"+str(i) for i in range(1, m+1)], 
					   columns=list(itertools.repeat("n_obs", ni.shape[1])))
		design["ni"] = ni


	### for model_switch ###
	if type(model_switch) is list:
		length = max([len(i) for i in model_switch])
		model_switch = np.array([np.pad(i, (0, length-len(i)), 'constant', constant_values=np.nan) for i in model_switch]) # convert a list of vectors to an array
	if model_switch is None:
		model_switch = xt * 0 + 1
	if size(model_switch)[0] == 1 and m != 1:
		model_switch  = np.tile(model_switch.flatten(), m).reshape(m, model_switch.size) # flatten xt, repeat by m times, and reshape to (col: xt's element number, row: m)

	if type(model_switch) is not np.ndarray and type(model_switch) is not pd.DataFrame:
		model_switch = np.array(model_switch)

	if test_mat_size(np.array(size(xt)), np.array(model_switch), "model_switch") == 1:
		model_switch = pd.DataFrame(model_switch, 
					   index=["grp_"+str(i) for i in range(1, m+1)], 
					   columns=["obs_"+str(i) for i in range(1, model_switch.shape[1]+1)])
		design["model_switch"] = model_switch


	### for a ###
	if a is not None:
		if type(a) == list:
			a = pd.DataFrame(a)
		elif len(a.shape) == 1:
			a = np.array([a])
		colnam = None
		if size(a)[0] == 1 and m != 1:
			a_ = []
			for i in range(0, a.shape[1]):
				for j in range(0, a.shape[0]):
					a_.append(a[j][i])
			a = pd.DataFrame(np.tile(a_, m).reshape(m, a.size),
										 columns=colnam, index=["grp_"+str(i) for i in range(1, m+1)])

		if type(a) is not np.ndarray and type(a) is not pd.DataFrame:
			a = np.array(a)

		a = pd.DataFrame(a)
		if size(a)[0] != m:
			raise Exception("The number of rows in a (" + str(size(a)[0]) + ") is not the same as the number of groups m (" + str(m) + ")")
		a.set_axis(["grp_"+str(i) for i in range(1, m+1)], axis="index")
		count = 0
		for i in range(0, a.shape[1]):
			if re.match(r'^X\d*$', str(a.columns[i])) != None:
				count += 1
		if count == size(a)[1]:
			a.set_axis([None] * a.shape[1], axis="column")
		design["a"] = a


	### for x ###
	if x != None:
		if type(x) == list:
			x = pd.DataFrame(x)
		colnam = x.columns.values.tolist()
		if size(x)[0] == 1 and m != 1:
			x_ = []
			for i in range(0, x.shape[1]):
				for j in range(0, x.shape[0]):
					x_.append(x[i][j])
			x = pd.DataFrame(np.tile(x_, m).reshape(m, x.size),
										 columns=colnam, index=["grp_"+str(i) for i in range(1, m+1)])

		if type(x) is not np.ndarray and type(x) is not pd.DataFrame:
			x = pd.DataFrame(x)

		if size(x)[0] != m:
			raise Exception("The number of rows in x (" + str(size(x)[0]) + "is not the same as the number of groups m (" + str(m) + ")")
		x.set_axis(["grp_"+str(i) for i in range(1, m+1)], axis="index")
		design["x"] = x


	### for groupsize ###
	if max(size(groupsize)) == 1 and m != 1:
		groupsize = [groupsize] * m
		groupsize = pd.DataFrame(np.array(groupsize).reshape([m ,1]), index=["grp_"+str(i) for i in range(1, m+1)])

		if type(groupsize) is not np.ndarray and type(groupsize) is not pd.DataFrame:
			groupsize = pd.DataFrame(groupsize)
			
		if test_mat_size(np.array([m, 1]), groupsize.to_numpy(), "groupsize") == 1:
			groupsize.set_axis(["grp_"+str(i) for i in range(1, m+1)], axis="index")
			groupsize.set_axis(["n_id"] * groupsize.shape[1], axis="columns")
		design["groupsize"] = groupsize
	

	return design