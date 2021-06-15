"""

Author: Caiya Zhang, Yuchen Zheng
"""


import project.all_modules as am


def create_design(
				  xt, # Matrix defining the initial sampling schedule row major, can also be a list of vectors
				  groupsize, # Vector defining the size of the different groups (num individuals in each group)
				  m = None, # Number of groups, computed from xt if not defined
				  x = None, # Matrix defining the initial discrete values
				  a = None,
				  ni = None, # Vector defining the number of samples for each group, computed as all elements of xt by default
				  model_switch = None # Vector defining which response a certain sampling time belongs to, defaults to one for all elements of xt
):

	design = {}

	### for xt, m ###
	if type(xt) is list: 
		length = max([len(i) for i in xt])
# int to float没写
		xt_ = []
		for i in range(0, len(xt)):
			xt[i] = xt[i].astype(am.np.float32)
			xt[i] = am.np.pad(xt[i], (0, length-len(xt[i])), "constant", constant_values=am.np.nan)
			# xt = am.np.array([am.np.pad(i, (0, length-len(i)), 'constant', constant_values=am.np.nan) for i in xt]) # convert a list of vectors to an array
			xt_.append(xt[i].tolist())
		xt = am.np.array(xt_)

	if m is None: 
		m = xt.shape[0] # get xt row (same as "m = am.size.size(xt, 1)")

	if am.size.size(xt, 1) == 1 and m != 1:
		xt = am.np.tile(xt.flatten(), m).reshape(m, xt.size) # flatten xt, repeat by m times, and reshape to (col: xt's element number, row: m)

# 没写！！！if(!is.matrix(xt)) xt <- rbind(xt)
	
	if (am.size.size(xt, 1) != m):
		raise Exception("The number of rows in xt (" + str(am.size.size(xt, 1)) + ") is not the same as the number of groups m (" + str(m) + ")")
	
	xt = am.pd.DataFrame(xt, 
					   index=["grp_"+str(i) for i in range(1, m+1)], 
					   columns=["obs_"+str(i) for i in range(1, xt.shape[1]+1)]) # same as "am.size.size(xt, 2)+1"
	
	design["xt"] = xt
# 没写！！！names(m) <- "n_grp"
	design["m"] = m


	### for ni ###
	if ni is None:
		ni = am.np.count_nonzero(1-am.np.isnan(xt), axis=1).reshape(xt.shape[0], 1)
	
	if type(ni) != am.np.ndarray:
		ni = am.np.array(ni).reshape([len(ni), 1])
	if am.test_mat_size.test_mat_size(am.np.array([m, 1]), ni, "ni") == 1:
		ni = am.pd.DataFrame(ni, 
					   index=["grp_"+str(i) for i in range(1, m+1)], 
					   columns=list(am.itertools.repeat("n_obs", ni.shape[1])))
		design["ni"] = ni


	### for model_switch ###
	if type(model_switch) is list:
		length = max([len(i) for i in model_switch])
		model_switch = am.np.array([am.np.pad(i, (0, length-len(i)), 'constant', constant_values=am.np.nan) for i in model_switch]) # convert a list of vectors to an array
	if model_switch is None:
		model_switch = xt * 0 + 1
	if am.size.size(model_switch, 1) == 1 and m != 1:
		model_switch  = am.np.tile(model_switch.flatten(), m).reshape(m, model_switch.size) # flatten xt, repeat by m times, and reshape to (col: xt's element number, row: m)
# 没写！！！if(!is.matrix(model_switch)) model_switch <- rbind(model_switch)
	if am.test_mat_size.test_mat_size(am.np.array(am.size.size(xt)), am.np.array(model_switch), "model_switch") == 1:
		model_switch = am.pd.DataFrame(model_switch, 
					   index=["grp_"+str(i) for i in range(1, m+1)], 
					   columns=["obs_"+str(i) for i in range(1, model_switch.shape[1]+1)])
		design["model_switch"] = model_switch


	### for a ###
	if a is not None:
		if type(a) == list:
			a = am.pd.DataFrame(a)
		elif len(a.shape) == 1:
			a = am.np.array([a])
		colnam = None
		if am.size.size(a, 1) == 1 and m != 1:
			a_ = []
			for i in range(0, a.shape[1]):
				for j in range(0, a.shape[0]):
					a_.append(a[j][i])
			a = am.pd.DataFrame(am.np.tile(a_, m).reshape(m, a.size),
										 columns=colnam, index=["grp_"+str(i) for i in range(1, m+1)])
# 没写！！！if(!is.matrix(a)) a <- rbind(a)
		a = am.pd.DataFrame(a)
		if am.size.size(a, 1) != m:
			raise Exception("The number of rows in a (" + str(am.size.size(a, 1)) + ") is not the same as the number of groups m (" + str(m) + ")")
		a.set_axis(["grp_"+str(i) for i in range(1, m+1)], axis="index")
		count = 0
		for i in range(0, a.shape[1]):
			if am.re.match(r'^X\d*$', str(a.columns[i])) != None:
				count += 1
		if count == am.size.size(a, 2):
			a.set_axis([None] * a.shape[1], axis="column")
		design["a"] = a


	### for x ###
	if x != None:
		if type(x) == list:
			x = am.pd.DataFrame(x)
		colnam = x.columns.values.tolist()
		if am.size.size(x, 1) == 1 and m != 1:
			x_ = []
			for i in range(0, x.shape[1]):
				for j in range(0, x.shape[0]):
					x_.append(x[i][j])
			x = am.pd.DataFrame(am.np.tile(x_, m).reshape(m, x.size),
										 columns=colnam, index=["grp_"+str(i) for i in range(1, m+1)])
# 没写！！！if(!is.matrix(x)) x <- rbind(x)
		if am.size.size(x, 1) != m:
			raise Exception("The number of rows in x (" + str(am.size.size(x, 1)) + "is not the same as the number of groups m (" + str(m) + ")")
		x.set_axis(["grp_"+str(i) for i in range(1, m+1)], axis="index")
		design["x"] = x


	### for groupsize ###
	if max(am.size.size(groupsize)) == 1 and m != 1:
		groupsize = [groupsize] * m
		groupsize = am.pd.DataFrame(am.np.array(groupsize).reshape([m ,1]), index=["grp_"+str(i) for i in range(1, m+1)])
# 没写！！！if(!is.matrix(groupsize)) groupsize <- rbind(groupsize)
		if am.test_mat_size.test_mat_size(am.np.array([m, 1]), groupsize.to_numpy(), "groupsize") == 1:
			groupsize.set_axis(["grp_"+str(i) for i in range(1, m+1)], axis="index")
			groupsize.set_axis(["n_id"] * groupsize.shape[1], axis="columns")
		design["groupsize"] = groupsize
	

	return design


# xt1 = [am.np.array([1.0, 2.0, 3.0]), am.np.array([1.0, 2.0, 3.0, 4.0])]


# design_1 = create_design(xt=xt1, groupsize=20)

# print(design_1)