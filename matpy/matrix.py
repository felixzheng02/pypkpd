"""
matpy is a self-defined class that allows users to 
name rows, columns, and elements of matrices
Authors: Caiya Zhang, Yuchen Zheng
"""


from logging import raiseExceptions
from typing import List
from xmlrpc.client import boolean
import numpy as np
import pandas as pd
from project.poped_choose import poped_choose
from project.cell import cell
from project.util import default_if_none


class Matrix:

	def __init__(self, data, shape: list = None, datanam: list = None, axisnam: list = None):
		"""
		@parameters:
			data: int, float, [int], [float], np.ndarray
			datanam: [str]
			axisnam: [[rowname], [colname]]
		"""

		# data field:
		# data: np.ndarray of Matrix
		# shape: the complete shape of the Matrix, including those of sub-matrices, 1-d is automatically transferred to 2-d
		# size: the complete size of the Matrix that includes sub-matrices
		# datanam: list of data name
		# axisnam: list of axis name, stored as [[str]]

		if type(data) is int or type(data) is float or type(data) is np.int32 or type(data) is np.int64 or type(data) is np.float32 or type(data) is np.float64:
			data = np.array([data]).reshape([1, 1])			
		elif type(data) is list:
			# needs to fill empty places by np.nan
			recursively_fill_list(data, np.nan)
			data = np.array(data)

		self.shape = select(shape, data.shape)
		if len(self.shape) == 1:
			self.shape = [1, self.shape[0]]
		self.data = data.reshape(tuple(self.shape))
		self.size = self.get_data().size
		self.datanam = datanam
		if self.datanam is not None:
			recursively_fill_list(self.datanam, None)
			self.datanam = np.array(self.datanam).reshape(self.shape).tolist()
		self.axisnam = axisnam

	def get_shape(self, index: int = None):
		"""
		get the shape of the Matrix, return as [int]
		"""
		if len(self.shape) == 1:
			return [1, self.shape[0]]
		return list(self.shape)

	# def get_

	def get_size(self) -> int:
		"""
		get the size of the Matrix, return as int
		"""
		return self.size

	def get_datanam(self):
		"""
		get the data names of the Matrix, return as [str]
		"""
		return self.datanam

	def get_axisnam(self, index: int = None):
		"""
		get the axis names of the Matrix, return as [[str]], 
		each [str] represents names in one axis
		"""
		if index is None:
			return self.axisnam
		else:
			if self.axisnam is None:
				return None
			else:
				return self.axisnam[index]

	def get_one_data(self, name: str = None, index: list = None):
		"""
		provided with either name (as str) or index (as list) of the data, return one data
		"""
		if name is None and index is None:
			raise Exception("Please specify the name or the index of the data.")

		if name is not None: # search by name
			if self.get_datanam() is not None:
				tmp_datanam = np.array(self.get_datanam()).flatten().tolist()
				for index in range(0, len(tmp_datanam)):
					if name == tmp_datanam[index]:
						return self.get_data().flatten()[index]
			raise Exception("'%s' does not exists." % name)
		elif index is not None: # search by index
			return self.get_data()[tuple(index)]
	
	def get_data(self):
		"""
		return all the data as np.ndarray
		"""
		return self.data

	def get_partial_matrix(self, index: list, axisnam: bool = True):
		"""
		only works for 2-d matrix as for now
		index: [[row_start, row_end], [col_start, col_end]]
		"""
		row_start = default_if_none(index[0][0], 0)
		row_end = default_if_none(index[0][1], self.get_shape(0))
		col_start = default_if_none(index[1][0], 0)
		col_end = default_if_none(index[1][1], self.get_shape(1))
		if self.get_axisnam(0) is None:
			rownam = None	
		else:
			rownam = self.get_axisnam(0)[row_start:row_end]
		if self.get_axisnam(1) is None:
			colnam = None	
		else:
			colnam = self.get_axisnam(1)[col_start:col_end]
		return Matrix(self.get_data()[row_start:row_end, col_start:col_end],
						axisnam=[rownam, colnam])

	def set_data(self, data, shape: list = None, datanam: bool = False, axisnam: bool = False):
		"""
		set data of the Matrix, datanam, colnam, rownam are lost
		@parameters:
			data: int, float, [int], [float], np.ndarray, Matrix
			shape: list
		"""
		if type(data) is int or type(data) is float or type(data) is np.int32 or type(data) is np.int64 or type(data) is np.float32 or type(data) is np.float64:
			data = np.array([data]).reshape([1, 1])
		elif type(data) is list:
			data = np.array(data)
		elif type(data) is np.ndarray:
			data = data
		else:
			raise Exception("must provide int, float, list, or np.ndarray as data input")
		self.shape = select(shape, data.shape)
		self.data = data.reshape(self.get_shape())
		if not datanam:
			self.datanam = None
		if not axisnam:
			self.axisnam = None

	def set_shape(self, shape: list, axisnam: list = None):
		"""
		set shape of the Matrix, keep axisnam if applicable
		"""
		self.shape = shape
		self.data = self.get_data().reshape(shape)
		if axisnam is not None:
			self.set_axisnam(axisnam)
			if len(self.get_axisnam()) != len(self.get_shape()):
				self.set_axisnam(None)

	def set_datanam(self, datanam: list):
		self.datanam = datanam

	def set_axisnam(self, axisnam: list):
		self.axisnam = axisnam
	
	def set_one_axisnam(self, index: int, new_axisnam: list):
		self.axisnam[index] = new_axisnam
			
	def set_one_data(self, new_data, new_datanam = None, name: str = None, index: list = None): # could only set by int or float
		"""
		provided with either name (as str) or index (as list) of the data that needs to be changed, 
		change the data and its data name
		"""
		if name is None and index is None:
			raise Exception("Please specify the name or the index of the data that needs to be changed.")

		if name is not None: # search by name
			if self.get_datanam() is not None:
				data = self.get_data().flatten()
				datanam = np.array(self.get_datanam()).flatten().tolist()
				for index in range(0, len(datanam)):
					if name == datanam[index]:
						data[index] = new_data
						self.data = data.reshape(self.get_shape())
						self.set_one_datanam(new_datanam, name=name)
						# datanam[index] = new_datanam
						# self.datanam = np.array(datanam).reshape(self.get_shape()).tolist()
						return
			raise Exception("'%s' does not exists." % name)
		elif index is not None: # search by index
			self.data[tuple(index)] = new_data
			if new_datanam is not None:
				self.set_one_datanam(new_datanam, index=index)
	
	def set_one_datanam(self, new_datanam, name: str = None, index: list = None):
		"""
		provided with either name (as str) or index (as list) of the data whose its name needs to be changed,
		change the data name
		"""
		if name is None and index is None:
			raise Exception("Please specify the name or the index of the data name that needs to be changed.")

		if name is not None: # search by name
			if self.get_datanam() is not None:
				datanam = np.array(self.get_datanam()).flatten().tolist()
				for index in range(0, len(datanam)):
					if name == datanam[index]:
						datanam[index] = new_datanam
						self.datanam = np.array(datanam).reshape(self.get_shape()).tolist()
						return
			raise Exception("'%s' does not exists." % name)
		elif index is not None: # search by index
			datanam = np.array(self.get_datanam())
			datanam[tuple(index)] = new_datanam
			self.set_datanam(datanam.tolist())

	def get_dataframe(self):
		return pd.DataFrame(self.get_data(),
							columns=self.get_colnam(),
							index=self.get_rownam())


	def repeat(self, n: list[int], shape: list[int] = None, datanam: bool = False, axisnam: bool = False):
		"""
		repeat the Matrix along row or column, specified by the n list
		repeat n[0] along row and n[1] along column
		reshape as specified
		datanam: if datanam are repeated, all set to None if 0
		axisnam: if axisnam are repeated, all set to None if 0
		"""
		data = np.tile(self.get_data(), n)
		if shape is None:
			shape = data.shape
		if len(shape) == 1:
			shape = [1, shape[0]]
		data = data.reshape(tuple(shape))
		self.set_data(data, datanam=True, axisnam=True)
		self.set_shape(data.shape)
		self.size = data.size
		if datanam:
			if self.get_datanam() is not None:
				self.set_datanam(np.tile(self.get_datanam(), n).reshape(self.shape).tolist())
		else:
			self.datanam = None
		if axisnam:
			if self.get_axisnam() is not None:
				for i in range(0, len(self.get_axisnam())):
					if self.get_axisnam()[i] is not None:
						self.get_axisnam()[i] = self.get_axisnam()[i] * int(self.get_shape()[i]/shape[i])
		else:
			self.axisnam = None


	def expand(self, shape: list[int], fill_value = np.nan):
		"""
		fill the Matrix with fill_value to shape
		"""
		new_data = cell(shape, fill_value)
		new_data[0:self.get_shape()[0], 0:self.get_shape()[1]] = self.get_data()
		self.set_data(new_data)
		if self.get_axisnam() is not None:
			for i in range(0, len(shape)):
				self.set_one_axisnam(i, self.get_axisnam()[i].expand((shape[i]-self.get_shape()[i]) * None))
		if self.get_datanam() is not None:
			new_datanam = cell(shape, None)
			new_datanam[0:self.get_shape()[0], 0:self.get_shape()[1]] = self.get_datanam()
			self.set_datanam(new_datanam)
		self.set_shape(shape)


def select(input, default):
	"""
	if given input, set target variable as input, else, set target variable as default
	"""
	if input is None:
		return default
	else:
		return input


def recursively_fill_list(input_list, value):
	"""
	designed as private method
	recursively fill np.nan values to blank spaces in a list,
	which probably includes multiple layers of sublist
	"""
	if type(input_list[0]) is list:
		length = max([len(sub_list) for sub_list in input_list])
		for i in range(0, len(input_list)):
			recursively_fill_list(input_list[i], value)
			# fill
			input_list[i] = input_list[i] + [value] * (length - len(input_list[i]))

	else: # sub-list structure ends
		return
