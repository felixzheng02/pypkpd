"""

Author: Caiya Zhang, Yuchen Zheng
"""

def poped_choose(arg1, arg2, exp):
	try:
		arg1
	except NameError:
		if type(arg2) is str and exp == 1:
			raise ValueError(arg2)
		else:
			return arg2
	else:
		if arg1 is not None:
			return arg1
		elif type(arg2) is str and exp == 1:
			raise ValueError(arg2)
		else:
			return arg2