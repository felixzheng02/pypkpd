"""
matpy is a self-defined class that allows users to 
name rows, columns, and elements of matrices
Authors: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
import pandas as pd
from project.poped_choose import poped_choose

class matrix:

	def __init__(self, data, shape = None, datanam: list = None, colnam: list = None, rownam: list = None):
		# data field:
		# data
		# shape
		# size
		# datanam
		# colnam
		# rownam
		if type(data) is matrix:
			self.shape = poped_choose(shape, data.shape, 0)
			self.data = np.array(data.data).reshape(shape)
			self.size = data.size
			self.datanam = poped_choose(datanam, data.datanam, 0)
			self.colnam = poped_choose(colnam, data.colnam, 0)
			self.rownam = poped_choose(rownam, data.rownam, 0)

		else:
			if shape is None:
				self.shape = (1, data.size)
			else:
				self.shape = shape
			self.data = np.array(data).reshape(shape)
			self.size = self.get_data().size
			self.datanam = datanam
			self.colnam = colnam
			self.rownam = rownam
	
	def get_data(self):
		return self.data

	def get_shape(self):
		return list(self.shape)

	def get_size(self):
		return self.size

	def get_datanam(self):
		return self.datanam

	def get_colnam(self):
		return self.colnam

	def get_rownam(self):
		return self.rownam

	def get_one_data(self, name: str):
		if self.get_datanam() is not None:
			for index in range(0, len(self.get_datanam())):
				if name == self.get_datanam()[index]:
					return self.get_data().flatten()[index]
		raise Exception("'%s' does not exists." % name)

	def set_datanam(self, datanam: list):
		self.datanam = datanam

	def set_colnam(self, colnam: list):
		self.colnam = colnam

	def set_rownam(self, rownam: list):
		self.rownam = rownam

	def set_axisnam(self, colnam: list, rownam: list):
		self.colnam = colnam
		self.rownam = rownam

	def create_dataframe(self):
		return pd.DataFrame(self.get_data(),
							columns=self.get_colnam(),
							index=self.get_rownam())