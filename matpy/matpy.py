"""
matpy is a self-defined class that allows users to 
name rows, columns, and elements of matrices

Authors: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
import pandas as pd


class matpy:

	def __init__(self, data: np.ndarray, shape: tuple, datanam: list = None, colnam: list = None, rownam: list = None):
		self.data = data.reshape(shape)
		self.shape = shape
		self.datanam = datanam
		self.colnam = colnam
		self.rownam = rownam
	
	def get_data(self):
		return self.data

	def get_shape(self):
		return self.shape

	def get_datanam(self):
		return self.datanam

	def get_colnam(self):
		return self.colnam

	def get_rownam(self):
		return self.rownam

	def get_one_data(self, name: str):
		if self.datanam is not None:
			for index in range(0, len(self.datanam)):
				if name == self.datanam[index]:
					return self.data.flatten()[index]
		raise Exception("'%s' does not exists." % name)

	def set_datanam(self, datanam: list):
		self.datanam = datanam

	def set_colnam(self, colnam: list):
		self.colnam = colnam

	def set_rownam(self, rownam: list):
		self.rownam = rownam

	def create_dataframe(self):
		return pd.DataFrame(self.data,
							columns=self.colnam,
							index=self.rownam)