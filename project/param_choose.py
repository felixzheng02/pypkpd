"""
Get parameter from popedInput and replace it by None if key does not exist.

Author: Caiya Zhang, Yuchen Zheng
"""


def param_choose(input: dict, rep_val, exp: int, param: list, *argv):
	tmp_dict = input
	for k in argv:
		if k not in list(tmp_dict.keys()):
			if exp == 1:
				raise Exception(rep_val)
			else:
				return rep_val
		else:
			tmp_dict = tmp_dict[k]
	if argv[-1] not in param:
		return tmp_dict
