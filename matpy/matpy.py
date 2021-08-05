import numpy as np


class matpy:

	def __init__(self, data: np.ndarray, shape: tuple, datanam: list, colnam: list, rownam: list):
		self.data = data.reshape(shape)
		self.shape = shape
		self.datanam = datanam
		self.colnam = colnam
		self.rownam = rownam
	
	def get_data(self):
		return self.data

	def get_one_data(self, name: str):
		for index in range(0, len(self.datanam)):
			if name == self.datanam[index]:
				return self.data.flatten()[index]
