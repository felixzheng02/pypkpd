import path
import unittest
import numpy as np
from project.create_design import create_design
from project.create_design_space import create_design_space
from matpy.matrix import Matrix

class TestCreateDesign(unittest.TestCase):

    def test_create_design_space(self):

        design_1 = create_design(xt=Matrix([[1,2,3,4,5], [1,2,3,4]]),
                         groupsize=Matrix([50, 20]),
                         a=Matrix([[70, 1000],[35, 1000]], axisnam=[None, ["WT", "DOSE"]]))
        
        ############################################        
        
        ds_1 = create_design_space(design_1)

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_1_a = create_design_space(design_1, our_zero = 1e-5)

        self.assertTrue(np.array_equal(ds_1_a["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1_a["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1_a["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1_a["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1_a["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1_a["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1_a["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1_a["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1_a["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1_a["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1_a["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1_a["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1_a["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1_a["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1_a["design_space"]["use_grouped_xt"])

        ############################################

        ds_2 = create_design_space(design_1, maxni=10, maxxt=10, minxt=0)

        self.assertTrue(np.array_equal(ds_2["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5, 5, 5, 5, 5, 5], [1, 2, 3, 4, 5, 5, 5, 5, 5, 5]]), equal_nan=True))
        self.assertListEqual(ds_2["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7", "obs_8", "obs_9", "obs_10"]])

        self.assertEqual(ds_2["design"]["m"].get_value(), 2)
        self.assertEqual(ds_2["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_2["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_2["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_2["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]), equal_nan=True))
        self.assertListEqual(ds_2["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7", "obs_8", "obs_9", "obs_10"]])
        
        self.assertTrue(np.array_equal(ds_2["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_2["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_2["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_2["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_2["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_2["design_space"]["maxni"].get_data(), np.array([[10], [10]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["maxtotni"].get_data(), np.array([[20]]), equal_nan=True))
        self.assertIsNone(ds_2["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_2["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_2["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_2["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_2["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_2["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_2["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_2["design_space"]["maxxt"].get_data(), np.array([[10, 10, 10, 10, 10, 10, 10, 10, 10, 10], [10, 10, 10, 10, 10, 10, 10, 10, 10, 10]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7", "obs_8", "obs_9", "obs_10"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["minxt"].get_data(), np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7", "obs_8", "obs_9", "obs_10"]])

        self.assertTrue(np.array_equal(ds_2["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]]), equal_nan=True))
        self.assertListEqual(ds_2["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7", "obs_8", "obs_9", "obs_10"]])

        self.assertFalse(ds_2["design_space"]["use_grouped_xt"])

        ############################################

        ds_3 = create_design_space(design_1,maxni=10, mingroupsize=20, maxxt=10, minxt=0)

        self.assertTrue(np.array_equal(ds_3["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5, 5, 5, 5, 5, 5], [1, 2, 3, 4, 5, 5, 5, 5, 5, 5]]), equal_nan=True))
        self.assertListEqual(ds_3["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5", "obs_6", "obs_7", "obs_8", "obs_9", "obs_10"]])

        self.assertEqual(ds_3["design"]["m"].get_value(), 2)
        self.assertEqual(ds_3["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_3["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_3["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_3["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_3["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_3["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_3["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_3["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_3["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_3["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_3["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_3["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_3["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_3["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_3["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_3["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_3["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_3["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_3["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_3["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_3["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_3["design_space"]["use_grouped_xt"])

        ############################################

        ds_4 = create_design_space(design_1, maxa=Matrix([100,2000]))

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_5 = create_design_space(design_1, mina=Matrix([10,20]))

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        design_2 = create_design(xt=Matrix([[1,2,3,4,5], [1,2,3,4]]),
                         groupsize=Matrix([50, 20]),
                         a=Matrix([[70, 1000],[35, 1000]], axisnam=[None, ["WT", "DOSE"]]),
                         x=Matrix([[1, 100], [2, 200]], axisnam=[None, ["SEX", "DOSE_discrete"]]))

        ds_6 = create_design_space(design_2) 

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_7 = create_design_space(design_2,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]))

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_8 = create_design_space(design_2,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                grouped_xt=np.array([1, 2, 3, 4, 5]))

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_9 = create_design_space(design_2,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                use_grouped_xt=True)

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################
                               
        design_3 = create_design(xt=Matrix([[1,2,3,4,5], [1,2,3,4]]),
                         groupsize=Matrix([50, 20]),
                         a=Matrix([35, 1000], axisnam=[None, ["WT", "DOSE"]]),
                         x=Matrix([1, 100], axisnam=[None, ["SEX", "DOSE_discrete"]]))

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_10 = create_design_space(design_3,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                use_grouped_a=True)

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_11 = create_design_space(design_2,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                grouped_a=np.array([[1, 2], [3, 2]]))

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_12 = create_design_space(design_3,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                use_grouped_xt=True)

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_13 = create_design_space(design_3,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                grouped_x=np.array([1, 2], [3, 2]))

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        # seq_1 = 1:10
        seq_1 = list(range(1, 11))
        
        ds_14 = create_design_space(design_1,maxxt=10,minxt=0,
                                    xt_space = [seq_1,seq_1,seq_1,seq_1,seq_1])
        
        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        ds_15 = create_design_space(design_1,maxxt=10,minxt=0,xt_space = [seq_1])

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])

        ############################################

        # ds_16 = create_design_space(design_1,maxxt=10,minxt=0,xt_space = xt_space)

        ds_17 = create_design_space(design_1,a_space=[list(range(1, 100)), list(range(1000, 101000, 1000))])

        self.assertTrue(np.array_equal(ds_1["design"]["xt"].get_data(),
                        np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertEqual(ds_1["design"]["m"].get_value(), 2)
        self.assertEqual(ds_1["design"]["m"].get_name(), "n_grp")

        self.assertTrue(np.array_equal(ds_1["design"]["ni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["ni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["model_switch"].get_data(), np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["model_switch"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])
        
        self.assertTrue(np.array_equal(ds_1["design"]["a"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design"]["groupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design"]["groupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxa"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxa"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mina"].get_data(), np.array([[70, 1000], [35, 1000]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mina"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_a"].get_data(), np.array([[1, 2], [3, 4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_a"].get_axisnam(), [["grp_1", "grp_2"], ["WT", "DOSE"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_a"])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minni"].get_data(), np.array([[5], [4]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minni"].get_axisnam(), [["grp_1", "grp_2"], ["n_obs"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotni"].get_data(), np.array([[9]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotni"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxgroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxgroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["mingroupsize"].get_data(), np.array([[50], [20]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["mingroupsize"].get_axisnam(), [["grp_1", "grp_2"], ["n_id"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxtotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["maxtotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["mintotgroupsize"].get_data(), np.array([[70]]), equal_nan=True))
        self.assertIsNone(ds_1["design_space"]["mintotgroupsize"].get_axisnam())

        self.assertTrue(np.array_equal(ds_1["design_space"]["maxxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["maxxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["minxt"].get_data(), np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["minxt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertTrue(np.array_equal(ds_1["design_space"]["grouped_xt"].get_data(), np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, np.nan]]), equal_nan=True))
        self.assertListEqual(ds_1["design_space"]["grouped_xt"].get_axisnam(), [["grp_1", "grp_2"], ["obs_1", "obs_2", "obs_3", "obs_4", "obs_5"]])

        self.assertFalse(ds_1["design_space"]["use_grouped_xt"])


if __name__ == '__main__':
    unittest.main()