"""
Get parameter from popedInput and replace it by None if key does not exist.

Author: Caiya Zhang, Yuchen Zheng
"""


def param_choose(input: dict, rep_val, exp: int, param: list = None, kwargs: dict = None, *argv):
	
	# determine if parameter's value is given as function input
	if param is not None and kwargs is not None:
		if argv[-1] in param:
			return kwargs[argv[-1]]
	
	# look for keys in dicts
	tmp_dict = input
	for k in argv:
		if k not in list(tmp_dict.keys()):
			if exp == 1:
				raise Exception(rep_val)
			else:
				return rep_val
		else:
			tmp_dict = tmp_dict[k]
	return tmp_dict
