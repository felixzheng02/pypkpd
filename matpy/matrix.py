"""
matpy is a self-defined class that allows users to 
name rows, columns, and elements of matrices
Authors: Caiya Zhang, Yuchen Zheng
"""


from logging import raiseExceptions
from typing import List
from hamcrest import none
import numpy as np
import pandas as pd
from project.poped_choose import poped_choose


class Matrix:

	def __init__(self, data, shape: list = None, datanam: list = None, axisnam: list = None):
		"""
		@parameters:
			data: int, float, [int], [float], np.ndarray, Matrix
			shape: list
			datanam: [str]
			colnam: [str]
			rownam: [str]
		"""

		# data field:
		# data: np.ndarray of Matrix
		# shape: the complete shape of the Matrix, including those of sub-matrices, 1-d is automatically transferred to 2-d
		# size: the complete size of the Matrix that includes sub-matrices
		# datanam: list of data name
		# axisnam: list of axis name, stored as [[str]]

		if type(data) is Matrix: # axisnam is not kept
			self.shape = select(shape, data.get_shape())
			self.data = np.array(data.get_data()).reshape(shape)
			self.size = data.get_size()
			self.datanam = select(datanam, data.get_datanam())
			self.axisnam = axisnam

		else:
			if type(data) is int or type(data) is float or type(data) is np.int32 or type(data) is np.int64 or type(data) is np.float32 or type(data) is np.float64:
				data = np.array([data]).reshape([1, 1])
			elif type(data) is list:
				data = np.array(data)
			self.shape = None
			self.shape = select(shape, data.shape)
			self.data = np.array(data).reshape(self.shape)
			self.size = self.get_data().size
			self.datanam = datanam
			self.axisnam = axisnam

	def get_shape(self):
		if len(self.shape) == 1:
			return [1, self.shape[0]]
		return list(self.shape)

	def get_size(self):
		return self.size

	def get_datanam(self):
		return self.datanam

	def get_axisnam(self):
		return self.axisnam

	def get_one_data(self, name: str = None, index: list = None):
		if name is None and index is None:
			raise Exception("Please specify the name or the index of the data.")

		if name is not None: # search by name
			if self.get_datanam() is not None:
				for index in range(0, len(self.get_datanam())):
					if name == self.get_datanam()[index]:
						return self.get_data().flatten()[index]
			raise Exception("'%s' does not exists." % name)
		elif index is not None: # search by index
			return self.get_data()[tuple(index)]
	
	def get_data(self):
		return self.data

	def set_data(self, data, shape: list = None, datanam: list = None):
		"""
		set data of the matrix, datanam, colnam, rownam are lost
		@parameters:
			data: int, float, [int], [float], np.ndarray, Matrix
			shape: list
		"""
		if type(data) is Matrix:
			data = data.get_data()
		elif type(data) is int or type(data) is float or type(data) is np.int32 or type(data) is np.int64 or type(data) is np.float32 or type(data) is np.float64:
			data = np.array([data]).reshape([1, 1])
		elif type(data) is list:
			data = np.array(data)
		elif type(data) is np.ndarray:
			data = data
		else:
			raise Exception("must provide int, float, list, np.ndarray, or Matrix as data input")
		self.shape = select(shape, data.shape)
		self.data = data.reshape(self.get_shape())
		self.set_datanam(datanam)
		self.set_axisnam(None)

	def set_shape(self, shape: list, axisnam: list = None):
		"""
		set shape of the matrix, keep axisnam if applicable
		"""
		self.shape = shape
		self.data = self.get_data().reshape(shape)
		if axisnam is not None:
			self.set_axisnam(axisnam)
		elif len(self.get_axisnam()) != len(self.get_shape()):
			self.set_axisnam(None)

	def set_datanam(self, datanam: list):
		self.datanam = datanam

	def set_axisnam(self, axisnam: list):
		self.axisnam = axisnam

	def set_one_data(self, new_data, new_datanam, name: str = None, index: list = None): # could only set by int or float
		if name is None and index is None:
			raise Exception("Please specify the name or the index of the data that needs to be changed.")

		if name is not None: # search by name
			data = self.get_data().flatten()
			if self.get_datanam() is not None:
				for index in range(0, len(self.get_datanam())):
					if name == self.get_datanam()[index]:
						data[index] = new_data
						self.data = data.reshape(self.get_shape())
						if new_datanam is not None:
							self.datanam[index] = new_datanam
						return
			raise Exception("'%s' does not exists." % name)
		elif index is not None: # search by index
			self.data[tuple(index)] = new_data
			if new_datanam is not None:
				self.set_one_datanam(new_datanam, index=index)
	
	def set_one_datanam(self, new_datanam, name: str = None, index: list = None):
		if name is None and index is None:
			raise Exception("Please specify the name or the index of the data name that needs to be changed.")

		if name is not None: # search by name
			if self.get_datanam() is not None:
				for index in range(0, len(self.get_datanam())):
					if name == self.get_datanam()[index]:
						self.datanam[index] = new_datanam
						return
			raise Exception("'%s' does not exists." % name)
		elif index is not None: # search by index
			datanam_index = 0
			for i in range(0, len(index)):
				datanam_index += index[i] * int(np.prod(self.get_shape()[i+1:]))
			self.datanam[datanam_index] = new_datanam

	def get_dataframe(self):
		return pd.DataFrame(self.get_data(),
							columns=self.get_colnam(),
							index=self.get_rownam())

def select(input, default):
	# if given input, set target variable as input, else, set target variable as default
	if input is None:
		return default
	else:
		return input