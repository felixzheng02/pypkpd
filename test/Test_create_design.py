from tkinter.messagebox import NO
from venv import create
import path
import unittest
import numpy as np
from project.create_design import create_design
from matpy.matrix import Matrix


class TestCreateDesign(unittest.TestCase):

    def test_create_design(self):
        
        xt1 = Matrix([[1, 2, 3], [1, 2, 3, 4]])
        xt4 = Matrix([[1,2,3,4,5], [1,2,3,4]]) 
        xt2 = Matrix(np.array([[1,2,3,4], [1,2,3,4]]))
        xt3 = Matrix(np.array([1,2,3,4]))

        ################################################

        design_1 = create_design(xt=xt1, groupsize=20)

        self.assertTrue(np.array_equal(design_1["xt"].get_data(), np.array([[1, 2, 3, np.nan], [1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_1["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_1["m"].get_value(), 2)
        self.assertEqual(design_1["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_1["ni"].get_data(), np.array([[3], [4]]), equal_nan=True))
        self.assertListEqual(design_1["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_1["model_switch"].get_data(), np.array([[1, 1, 1, np.nan], [1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_1["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertTrue(np.array_equal(design_1["groupsize"].get_data(), np.array([[20], [20]]), equal_nan=True))
        self.assertListEqual(design_1["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
        
        ################################################

        design_2 = create_design(xt=xt4, groupsize=20)

        self.assertTrue(np.array_equal(design_2["xt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(design_2["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertEqual(design_2["m"].get_value(), 2)
        self.assertEqual(design_2["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_2["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(design_2["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_2["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(design_2["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(design_2["groupsize"].get_data(), np.array([[20], [20]]), equal_nan=True))
        self.assertListEqual(design_2["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################

        design_3 = create_design(xt=xt2, groupsize=20)

        self.assertTrue(np.array_equal(design_3["xt"].get_data(), np.array([[1, 2, 3, 4], [1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_3["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_3["m"].get_value(), 2)
        self.assertEqual(design_3["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_3["ni"].get_data(), np.array([[4], [4]]), equal_nan=True))
        self.assertListEqual(design_3["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_3["model_switch"].get_data(), np.array([[1, 1, 1, 1], [1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_3["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertTrue(np.array_equal(design_3["groupsize"].get_data(), np.array([[20], [20]]), equal_nan=True))
        self.assertListEqual(design_3["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################

        #design_4 duplicates
        design_4 = create_design(xt=xt3, groupsize=20)

        self.assertTrue(np.array_equal(design_4["xt"].get_data(), np.array([[1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_4["xt"].get_axisnam(), [["grp_1"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_4["m"].get_value(), 1)
        self.assertEqual(design_4["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_4["ni"].get_data(), np.array([[4]]), equal_nan=True))
        self.assertListEqual(design_4["ni"].get_axisnam(), [["grp_1"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_4["model_switch"].get_data(), np.array([[1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_4["model_switch"].get_axisnam(), [["grp_1"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertTrue(np.array_equal(design_4["groupsize"].get_data(), np.array([[20]]), equal_nan=True))
        self.assertListEqual(design_4["groupsize"].get_axisnam(), [["grp_1"], ["n_id"]])
       
        ################################################
        
        design_5 = create_design(xt=xt3, groupsize=20,m=3)
        
        self.assertTrue(np.array_equal(design_5["xt"].get_data(), np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_5["xt"].get_axisnam(), [["grp_1", "grp_2", "grp_3"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_5["m"].get_value(), 3)
        self.assertEqual(design_5["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_5["ni"].get_data(), np.array([[4], [4], [4]]), equal_nan=True))
        self.assertListEqual(design_5["ni"].get_axisnam(), [["grp_1", "grp_2", "grp_3"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_5["model_switch"].get_data(), np.array([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_5["model_switch"].get_axisnam(), [["grp_1", "grp_2", "grp_3"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertTrue(np.array_equal(design_5["groupsize"].get_data(), np.array([[20], [20], [20]]), equal_nan=True))
        self.assertListEqual(design_5["groupsize"].get_axisnam(), [["grp_1", "grp_2", "grp_3"], ["n_id"]])
       
        ################################################
        
        design_6 = create_design(xt=xt1, groupsize=20, model_switch=Matrix(np.ones([2, 4])))
        
        self.assertTrue(np.array_equal(design_6["xt"].get_data(), np.array([[1, 2, 3, np.nan], [1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_6["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_6["m"].get_value(), 2)
        self.assertEqual(design_6["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_6["ni"].get_data(), np.array([[3], [4]]), equal_nan=True))
        self.assertListEqual(design_6["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_6["model_switch"].get_data(), np.array([[1, 1, 1, 1], [1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_6["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertTrue(np.array_equal(design_6["groupsize"].get_data(), np.array([[20], [20]]), equal_nan=True))
        self.assertListEqual(design_6["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################
        
        design_7 = create_design(xt=xt1, groupsize=20, a=Matrix(np.array([2, 3, 4])))
        
        self.assertTrue(np.array_equal(design_7["xt"].get_data(), np.array([[1, 2, 3, np.nan], [1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_7["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_7["m"].get_value(), 2)
        self.assertEqual(design_7["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_7["ni"].get_data(), np.array([[3], [4]]), equal_nan=True))
        self.assertListEqual(design_7["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_7["model_switch"].get_data(), np.array([[1, 1, 1, np.nan], [1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_7["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])

        self.assertTrue(np.array_equal(design_7["a"].get_data(), np.array([[2, 3, 4], [2, 3, 4]]), equal_nan=True))
        #self.assertListEqual(design_7["a"].get_axisnam(), [["grp_1", "grp_2"], [",1", ",2", ",3"]])

        self.assertTrue(np.array_equal(design_7["groupsize"].get_data(), np.array([[20], [20]]), equal_nan=True))
        self.assertListEqual(design_7["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################
        
        design_8 = create_design(xt=xt1, groupsize=20, a=Matrix(np.array([[2, 3, 4], [4, 5, 6]])))

        self.assertTrue(np.array_equal(design_8["xt"].get_data(), np.array([[1, 2, 3, np.nan], [1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_8["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_8["m"].get_value(), 2)
        self.assertEqual(design_8["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_8["ni"].get_data(), np.array([[3], [4]]), equal_nan=True))
        self.assertListEqual(design_8["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_8["model_switch"].get_data(), np.array([[1, 1, 1, np.nan], [1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_8["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertTrue(np.array_equal(design_8["a"].get_data(), np.array([[2, 3, 4], [4, 5, 6]]), equal_nan=True))
        #self.assertListEqual(design_8["a"].get_axisnam(), [["grp_1", "grp_2"], [",1", ",2", ",3"]])

        self.assertTrue(np.array_equal(design_8["groupsize"].get_data(), np.array([[20], [20]]), equal_nan=True))
        self.assertListEqual(design_8["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################
        
        # design_9 warnings and "a" duplicates
        design_9 = create_design(xt=xt1, groupsize=20, a=[[2, 3, 4, 6], [4, 5, 6]])
        
        self.assertTrue(np.array_equal(design_9["xt"].get_data(), np.array([[1, 2, 3, np.nan], [1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_9["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_9["m"].get_value(), 2)
        self.assertEqual(design_9["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_9["ni"].get_data(), np.array([[3], [4]]), equal_nan=True))
        self.assertListEqual(design_9["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_9["model_switch"].get_data(), np.array([[1, 1, 1, np.nan], [1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_9["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertTrue(np.array_equal(design_9["a"].get_data(), np.array([[2, 3, 4, 6], [4, 5, 6, np.nan]]), equal_nan=True))
        #self.assertListEqual(design_9["a"].get_axisnam(), [["grp_1", "grp_2"], [",1", ",2", ",3", ",4"]])

        self.assertTrue(np.array_equal(design_9["groupsize"].get_data(), np.array([[20], [20]]), equal_nan=True))
        self.assertListEqual(design_9["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################
        
        design_10 = create_design(xt=xt1, groupsize=20, a=[[2, 3, 4], [4, 5, 6]])

        self.assertTrue(np.array_equal(design_10["xt"].get_data(), np.array([[1, 2, 3, np.nan], [1, 2, 3, 4]]), equal_nan=True))
        self.assertListEqual(design_10["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertEqual(design_10["m"].get_value(), 2)
        self.assertEqual(design_10["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_10["ni"].get_data(), np.array([[3], [4]]), equal_nan=True))
        self.assertListEqual(design_10["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_10["model_switch"].get_data(), np.array([[1, 1, 1, np.nan], [1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_10["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4"]])
        
        self.assertTrue(np.array_equal(design_10["a"].get_data(), np.array([[2, 3, 4], [4, 5, 6]]), equal_nan=True))
        #self.assertListEqual(design_10["a"].get_axisnam(), [["grp_1", "grp_2"], [",1", ",2", ",3"]])         

        self.assertTrue(np.array_equal(design_10["groupsize"].get_data(), np.array([[20], [20]]), equal_nan=True))
        self.assertListEqual(design_10["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################
        
        design_11 = create_design(xt=Matrix([0, 1, 2, 4, 6, 8, 24]), groupsize=50, 
                                    a=Matrix([70, 1000], axisnam=[None, ["WT", "DOSE"]]))
        
        self.assertTrue(np.array_equal(design_11["xt"].get_data(), np.array([[0, 1, 2, 4, 6, 8, 24]]), equal_nan=True))
        self.assertListEqual(design_11["xt"].get_axisnam(), [["grp_1"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7"]])
        
        self.assertEqual(design_11["m"].get_value(), 1)
        self.assertEqual(design_11["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_11["ni"].get_data(), np.array([[7]]), equal_nan=True))
        self.assertListEqual(design_11["ni"].get_axisnam(), [["grp_1"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_11["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_11["model_switch"].get_axisnam(), [["grp_1"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7"]])
        
        self.assertTrue(np.array_equal(design_11["a"].get_data(), np.array([[70, 1000]]), equal_nan=True))
        self.assertListEqual(design_11["a"].get_axisnam(), [["grp_1"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(design_11["groupsize"].get_data(), np.array([[50]]), equal_nan=True))
        self.assertListEqual(design_11["groupsize"].get_axisnam(), [["grp_1"], ["n_id"]])
       
        ################################################
        
        design_12 = create_design(xt=Matrix([0, 1, 2, 4, 6, 8, 24]), groupsize=50, 
                                    a=Matrix([70, 1000], axisnam=[None, ["WT", "DOSE"]]), m=2)

        self.assertTrue(np.array_equal(design_12["xt"].get_data(), np.array([[0, 1, 2, 4, 6, 8, 24], [0, 1, 2, 4, 6, 8, 24]]), equal_nan=True))
        self.assertListEqual(design_12["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7"]])
        
        self.assertEqual(design_12["m"].get_value(), 2)
        self.assertEqual(design_12["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_12["ni"].get_data(), np.array([[7], [7]]), equal_nan=True))
        self.assertListEqual(design_12["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_12["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_12["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7"]])
        
        self.assertTrue(np.array_equal(design_12["a"].get_data(), np.array([[70, 1000], [70, 1000]]), equal_nan=True))
        self.assertListEqual(design_12["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(design_12["groupsize"].get_data(), np.array([[50], [50]]), equal_nan=True))
        self.assertListEqual(design_12["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################
        
        design_13 = create_design(xt=Matrix([0, 1, 2, 4, 6, 8, 24]), groupsize=50, 
                                    a=Matrix([[70, 1000], [200, 90, 45]], axisnam=[None, ["WT", "DOSE", "AGE"]]), m=2)

        self.assertTrue(np.array_equal(design_13["xt"].get_data(), np.array([[0, 1, 2, 4, 6, 8, 24], [0, 1, 2, 4, 6, 8, 24]]), equal_nan=True))
        self.assertListEqual(design_13["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7"]])
        
        self.assertEqual(design_13["m"].get_value(), 2)
        self.assertEqual(design_13["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_13["ni"].get_data(), np.array([[7], [7]]), equal_nan=True))
        self.assertListEqual(design_13["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_13["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(design_13["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7"]])
        
        self.assertTrue(np.array_equal(design_13["a"].get_data(), np.array([[70, 1000, np.nan], [200, 90, 45]]), equal_nan=True))
        self.assertListEqual(design_13["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE", "AGE"]])

        self.assertTrue(np.array_equal(design_13["groupsize"].get_data(), np.array([[50], [50]]), equal_nan=True))
        self.assertListEqual(design_13["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################
        
        # design_14 same as design_13

        design_15 = create_design(xt=xt4, groupsize=np.array([50, 20]), 
                                    a=Matrix([[2, 3, 4], [4, 5, 6]], axisnam=[None, ["DOSE", "WT", "AGE"]]))
        
        self.assertTrue(np.array_equal(design_15["xt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(design_15["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertEqual(design_15["m"].get_value(), 2)
        self.assertEqual(design_15["m"].get_name(), "n_grp")
        
        self.assertTrue(np.array_equal(design_15["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(design_15["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(design_15["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(design_15["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(design_15["a"].get_data(), np.array([[2, 3, 4], [4, 5, 6]]), equal_nan=True))
        self.assertListEqual(design_15["a"].get_axisnam(), [["grp_1", "grp_2"], ["DOSE", "WT", "AGE"]])

        self.assertTrue(np.array_equal(design_15["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(design_15["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])
       
        ################################################
        
if __name__ == '__main__':
    unittest.main()