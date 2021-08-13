import numpy as np
from matpy.matrix import matrix


data_1 = matrix(data=[matrix(np.array([1, 2, 3])), matrix(np.array([1, 2, 4]))])
print(data_1.get_all_data().reshape(2, 3))

# data_2 = matrix(np.array([1, 2, 3, 4]),
# 			 (2, 2),
# 			 ["one", "two", "three", "four"],
# 			 ["col_1", "col_2"],
# 			 ["row_1", "row_2"])
# data_2.set_one_data(matrix(np.array([1, 2, 3])), index=[0, 0])
# print(data_2.get_data()[0][0].get_data())

# print(data.get_data())
# print(data.get_one_data("one"))
# print(matrix(data).get_data())
# data.set_multiple_data(0, [0, 1], [0, 1])
# print(data.get_data())
# data.set_one_data(np.array([1, 2, 3]), index=[0, 0])