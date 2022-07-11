"""
For parameter choose in create_poped_database()

Author: Caiya Zhang, Yuchen Zheng
"""


def param_choose(param, pypkpdInput, replace, key_1, key_2=None, key_3=None):
	if param is not None:
		return param
	if pypkpdInput is not None:
		if key_1 is not None:
			if key_1 in pypkpdInput.keys():
				if key_2 is not None:
					if key_2 in pypkpdInput[key_1].keys():
						if key_3 is not None:
							if key_3 in pypkpdInput[key_1][key_2].keys():
								return pypkpdInput[key_1][key_2][key_3]
							else: raise Exception(key_3 + "not in keys.")
						else: return pypkpdInput[key_1][key_2]
					else: raise Exception(key_2 + "not in keys.")
				else: return pypkpdInput[key_1]
			else: raise Exception(key_1 + "not in keys.")
		else:
			return Exception("Must provide at least one key.")
	raise Exception("Must provide pypkpdInput dict.")
