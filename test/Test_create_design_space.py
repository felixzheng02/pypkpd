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
        
        ds_1 = create_design_space(design_1)

        ds_1_a = create_design_space(design_1, our_zero = 1e-5)

        ds_2 = create_design_space(design_1, maxni=10, maxxt=10, minxt=0)

        ds_3 = create_design_space(design_1,maxni=10, mingroupsize=20, maxxt=10, minxt=0)

        ds_4 = create_design_space(design_1, maxa=Matrix([100,2000]))

        ds_5 = create_design_space(design_1, mina=Matrix([10,20]))

        design_2 = create_design(xt=Matrix([[1,2,3,4,5], [1,2,3,4]]),
                         groupsize=Matrix([50, 20]),
                         a=Matrix([[70, 1000],[35, 1000]], axisnam=[None, ["WT", "DOSE"]]),
                         x=Matrix([[1, 100], [2, 200]], axisnam=[None, ["SEX", "DOSE_discrete"]]))

        ds_6 = create_design_space(design_2) 

        ds_7 = create_design_space(design_2,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]))

        ds_8 = create_design_space(design_2,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                grouped_xt=np.array([1, 2, 3, 4, 5]))

        ds_9 = create_design_space(design_2,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                use_grouped_xt=True)
                               
        design_3 = create_design(xt=Matrix([[1,2,3,4,5], [1,2,3,4]]),
                         groupsize=Matrix([50, 20]),
                         a=Matrix([35, 1000], axisnam=[None, ["WT", "DOSE"]]),
                         x=Matrix([1, 100], axisnam=[None, ["SEX", "DOSE_discrete"]]))

        ds_10 = create_design_space(design_3,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                use_grouped_a=True)

        ds_11 = create_design_space(design_2,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                grouped_a=np.array([[1, 2], [3, 2]]))

        ds_12 = create_design_space(design_3,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                use_grouped_xt=True)

        ds_13 = create_design_space(design_3,
                                x_space=Matrix([[1, 2], list(range(100, 420, 20))], 
                                axisnam=[["SEX", "DOSE_discrete"], None]),
                                grouped_x=np.array([1, 2], [3, 2]))

        # seq_1 = 1:10
        seq_1 = list(range(1, 11))
        
        ds_14 = create_design_space(design_1,maxxt=10,minxt=0,
                                    xt_space = [seq_1,seq_1,seq_1,seq_1,seq_1])
        ds_15 = create_design_space(design_1,maxxt=10,minxt=0,xt_space = [seq_1])

        # ds_16 = create_design_space(design_1,maxxt=10,minxt=0,xt_space = xt_space)

        ds_17 = create_design_space(design_1,a_space=[list(range(1, 100)), list(range(1000, 101000, 1000))])

if __name__ == '__main__':
    unittest.main()