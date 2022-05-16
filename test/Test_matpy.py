import unittest
import numpy as np
from matpy.matrix import Matrix


class TestMatpy(unittest.TestCase):

    def test_matrix_init(self):
        
        matrix_1 = Matrix(
                        np.array([1, 2, 3, 4]),
			            [2, 2],
			            ["one", "two", "three", "four"],
			            [["row_1", "row_2"], ["col_1", "col_2"]]
                        )
        matrix_2 = Matrix(1)
        matrix_3 = Matrix(2.0)
        matrix_4 = Matrix([10, 20, 30.0, 40], [1, 1, 4])

        self.assertListEqual(matrix_1.get_shape(), [2, 2])
        self.assertListEqual(matrix_2.get_shape(), [1, 1])
        self.assertListEqual(matrix_3.get_shape(), [1, 1])
        self.assertListEqual(matrix_4.get_shape(), [1, 1, 4])
        print("get_shape() test passed")
        self.assertEqual(matrix_1.get_size(), 4)
        self.assertEqual(matrix_2.get_size(), 1)
        self.assertEqual(matrix_3.get_size(), 1)
        self.assertEqual(matrix_4.get_size(), 4)
        print("get_size() test passed")
        self.assertListEqual(matrix_1.get_datanam(), ["one", "two", "three", "four"])
        print("get_datanam() test passed")
        self.assertListEqual(matrix_1.get_axisnam()[0], ["row_1", "row_2"])
        print("get_axisnam() test passed")
        self.assertEqual(matrix_1.get_one_data(name="two"), 2)
        self.assertEqual(matrix_1.get_one_data(index=[1, 1]), 4)
        self.assertEqual(matrix_2.get_one_data(index=[0, 0]), 1)
        self.assertEqual(matrix_3.get_one_data(index=[0, 0]), 2)
        self.assertEqual(matrix_4.get_one_data(index=[0, 0, 3]), 40.0)
        print("get_one_data() test passed")
        print("test 1 passed")

    def test_set_data(self):
        matrix = Matrix(
                        np.array([1, 2, 3, 4]),
			            [2, 2],
			            ["one", "two", "three", "four"],
			            [["row_1", "row_2"], ["col_1", "col_2"]]
                        )
        matrix.set_data(1)
        self.assertListEqual(matrix.get_shape(), [1, 1])
        matrix.set_data(np.array([1, 2, 3, 4]), [1, 2, 2])
        self.assertEqual(matrix.get_one_data(index=[0, 0, 1]), 2)
        matrix.set_data([10, 20, 30], [3, 1])
        self.assertEqual(matrix.get_one_data(index=[2, 0]), 30)
        matrix.set_data(matrix, [1, 3])
        self.assertListEqual(matrix.get_shape(), [1, 3])
        print("test 2 passed")
    
    def test_set_shape(self):
        matrix = Matrix(
                        np.array([1, 2, 3, 4]),
			            [2, 2],
			            ["one", "two", "three", "four"],
			            [["row_1", "row_2"], ["col_1", "col_2"]]
                        )
        matrix.set_shape([1, 2, 2], [["axis1"], ["axis21", "axis22"], ["axis31", "axis32"]])
        self.assertListEqual(matrix.get_shape(), [1, 2, 2])
        self.assertListEqual(matrix.get_axisnam()[0], ["axis1"])
        self.assertEqual(matrix.get_one_data(index=[0, 1, 1]), 4)
        self.assertEqual(matrix.get_one_data(name="one"), 1)
        print("test 3 passed")
    
    def test_set_datanam(self):
        matrix = Matrix(
                        np.array([1, 2, 3, 4]),
			            [2, 2],
			            ["one", "two", "three", "four"],
			            [["row_1", "row_2"], ["col_1", "col_2"]]
                        )
        matrix.set_shape([1, 4, 1])
        matrix.set_datanam(["1", "2", "3", "4"])
        self.assertEqual(matrix.get_one_data(name="1"), 1)
        print("test 4 passed")
    
    def test_set_axisnam(self):
        matrix = Matrix(
                        np.array([1, 2, 3, 4]),
			            [2, 2],
			            ["one", "two", "three", "four"],
			            [["row_1", "row_2"], ["col_1", "col_2"]]
                        )
        matrix.set_axisnam([["axis11", "axis12"], ["axis21", "axis22"]])
        self.assertListEqual(matrix.get_axisnam()[1], ["axis21", "axis22"])
        print("test 5 passed")

    def test_set_one_data(self):
        matrix = Matrix(
                        np.array([1, 2, 3, 4]),
			            [2, 2],
			            ["one", "two", "three", "four"],
			            [["row_1", "row_2"], ["col_1", "col_2"]]
                        )
        matrix.set_one_data(20, "twenty", name="two")
        self.assertEqual(matrix.get_one_data(index=[0, 1]), 20)
        matrix.set_one_data(30, "thirty", index=[1, 0])
        matrix.set_shape([1, 2, 2])
        self.assertEqual(matrix.get_one_data(name="thirty"), 30)
        self.assertEqual(matrix.get_one_data(index=[0, 1, 0]), 30)
        print("test 6 passed")

if __name__ == '__main__':
    unittest.main()