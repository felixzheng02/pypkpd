import numpy as np


class num:

	def __init__(self, input):
		# data field:
		# value
		# type
		self.value = input
		self.type = type(self.value)
	
	def get_value(self):
		return self.value

	def get_type(self):
		return self.type

	def trans_array(self):
		"""return self.value as an np.ndarray"""
		if self.get_type() is int or self.get_type() is float:
			return np.array([self.get_value()])
		else: # self.get_type() is list or np.ndarray
			return np.array(self.get_value())

	def get_by_index(self, index: int): # only 1-d array is considered
		"""return value by index, nan if out of index"""
		index_max = np.size(self.trans_array())
		if index < index_max:
			return self.trans_array()[index]
		else:
			return np.nan
