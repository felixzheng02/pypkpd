"""
matpy is a self-defined class that allows users to 
name rows, columns, and elements of matrices
Authors: Caiya Zhang, Yuchen Zheng
"""


from typing import List
import numpy as np
import pandas as pd
from project.poped_choose import poped_choose


class matrix:

	# shape parameter does not function well
	def __init__(self, data, shape = None, datanam: list = None, colnam: list = None, rownam: list = None):
		# data field:
		# data: only represents the top-level structure of the matrix
		# shape: the complete shape of the matrix, including those of sub-matrices
		# size: the complete size of the matrix that includes sub-matrices
		# datanam
		# colnam
		# rownam
		if type(data) is matrix:
			self.shape = poped_choose(shape, data.get_shape(), 0)
			self.data = np.array(data.get_data()).reshape(shape)
			self.size = data.get_size()
			self.datanam = poped_choose(datanam, data.datanam, 0)
			self.colnam = poped_choose(colnam, data.colnam, 0)
			self.rownam = poped_choose(rownam, data.rownam, 0)

		elif type(data) is int or type(data) is float or type(data) is np.int64 or type(data) is np.float64:
			self.shape = [1, 1]
			self.data = np.array([data])
			self.size = 1
			self.datanam = datanam
			self.colnam = colnam
			self.rownam = rownam

		elif all(isinstance(n, matrix) for n in data) and type(data) is list:
			self.data = np.array(data)
			tmp_shape = [len(data)]
			for i in self.get_data()[0].get_shape():
				if i != 1:
					tmp_shape.append(i)
			self.shape = poped_choose(shape, tmp_shape, 0)
			self.size = len(data)
			self.datanam = datanam
			self.colnam = colnam
			self.rownam = rownam

		elif type(data) is np.ndarray:
			self.shape = poped_choose(shape, data.shape, 0)
			self.data = np.array(data).reshape(shape)
			# if len(self.get_shape()) == 1:
			# 	self.shape = [1, self.get_shape()[0]]
			# if shape is None:
			# 	self.shape = data.shape
			# 	if len(self.shape) == 1:
			# 		self.shape = [1, self.shape[0]]
			# else:
			# 	self.shape = shape
			self.size = self.get_all_data().size
			self.datanam = datanam
			self.colnam = colnam
			self.rownam = rownam
	
	def get_data(self):
		return self.data

	def get_shape(self):
		if len(self.shape) == 1:
			return [1, self.shape[0]]
		return list(self.shape)

	def get_size(self):
		return self.size

	def get_datanam(self):
		return self.datanam

	def get_colnam(self):
		return self.colnam

	def get_rownam(self):
		return self.rownam

	def get_one_data(self, name: str = None, index = None):
		if name is not None:
			if self.get_datanam() is not None:
				for index in range(0, len(self.get_datanam())):
					if name == self.get_datanam()[index]:
						return self.get_all_data().flatten()[index]
			raise Exception("'%s' does not exists." % name)
		elif index is not None:
			return self.get_all_data()[list(index)[0]][list(index)[1]]
		else:
			raise Exception("Please specify the name or the index of the data needed.")
	
	def get_all_data(self):
		tmp = self
		if tmp.get_data().size == 0:
			return np.array([])
		elif type(tmp.get_data()[0]) is matrix:
			length = len(tmp.get_data().tolist())
			result = [None] * length
			for index in range(0, length):
				result[index] = tmp.get_data().tolist()[index].get_all_data()
			return np.array(result)
		else:
			return tmp.get_data()

	def set_data(self, data):
		self.data = np.array(data).reshape(self.get_shape())

	def set_shape(self, shape, colnam: None, rownam: None):
		self.shape = shape
		self.data = self.get_data().reshape(shape)
		if colnam is not None:
			self.set_colnam(colnam)
		if rownam is not None:
			self.set_rownam(rownam)
		if self.get_colnam() is not None and self.get_colnam() is not None:
			if (len(self.get_rownam()) != len(self.get_shape()[0]) or 
			len(self.get_colnam()) != len(self.get_shape()[1])):
				raise Exception("Number of column names or row names does not match the given shape.")

	def set_datanam(self, datanam: list):
		self.datanam = datanam

	def set_colnam(self, colnam: list):
		self.colnam = colnam

	def set_rownam(self, rownam: list):
		self.rownam = rownam

	def set_axisnam(self, colnam: list, rownam: list):
		self.colnam = colnam
		self.rownam = rownam
	
	def set_one_data(self, new_data, name: str = None, index = None):
		if name is not None:
			data = self.get_data().reshape([self.get_data().size, ]).tolist()
			if self.get_datanam() is not None:
				for index in range(0, len(self.get_datanam())):
					if name == self.get_datanam()[index]:
						data[index] = new_data
						self.data = np.array(data).reshape(self.get_shape())
			raise Exception("'%s' does not exists." % name)
		elif index is not None:
			if type(new_data) is matrix:
				data = self.get_data().reshape([self.get_data().size, ]).tolist()
				data[(list(index)[0]+1) * (list(index)[1]+1) - 1] = new_data
				self.data = np.array(data).reshape(self.get_shape())
			else:
				self.get_data()[list(index)[0]][list(index)[1]] = new_data
		else:
			raise Exception("Please specify the name or the index of the data that needs to be changed.")

	def set_multiple_data(self, new_data, row = None, col = None): # row, col: [init, end]
		if row is None:
			self.get_data()[:, list(col)[0]:list(col)[1]] = new_data
		elif col is None:
			self.get_data()[list(row)[0]:list(row)[1], :] = new_data
		elif row is not None and col is not None:
			self.get_data()[list(row)[0]:list(row)[1], list(col)[0]:list(col)[1]] = new_data

	def create_dataframe(self):
		return pd.DataFrame(self.get_data(),
							columns=self.get_colnam(),
							index=self.get_rownam())