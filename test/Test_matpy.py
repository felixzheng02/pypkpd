import numpy as np
from matpy.matpy import matpy


data = matpy(np.array([1, 2, 3, 4]),
			 (2, 2),
			 ["one", "two", "three", "four"],
			 ["col_1", "col_2"],
			 ["row_1", "row_2"])

print(data.get_data())
print(data.get_one_data("one"))