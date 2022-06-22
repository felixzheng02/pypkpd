from tkinter.messagebox import NO
from venv import create
import path
import unittest
import numpy as np
from project.create_design import create_design
from project.ones import ones
from matpy.matrix import Matrix


class TestCreateDesign(unittest.TestCase):

    def test_create_design_init(self):
        
        xt1 = Matrix([[1, 2, 3], [1, 2, 3, 4]])
        xt4 = Matrix([[1,2,3,4,5], [1,2,3,4]]) 
        xt2 = Matrix(np.array([[1,2,3,4], [1,2,3,4]]))
        xt3 = Matrix(np.array([1,2,3,4]))

        design_1 = create_design(xt=xt1, groupsize=20)
        design_2 = create_design(xt=xt4, groupsize=20)
        design_3 = create_design(xt=xt2, groupsize=20)

        #design_4 duplicates
        design_4 = create_design(xt=xt3, groupsize=20)

        design_5 = create_design(xt=xt3, groupsize=20,m=3)
        design_6 = create_design(xt=xt1, groupsize=20, model_switch=ones([2, 4]))
        design_7 = create_design(xt=xt1, groupsize=20, a=Matrix(np.array([2, 3, 4])))
        design_8 = create_design(xt=xt1, groupsize=20, a=Matrix(np.array([[2,3,4], [4,5,6]])))

        # design_9 warnings and "a" duplicates
        design_9 = create_design(xt=xt1, groupsize=20, a=[[2,3,4,6], [4,5,6]])
        design_10 = create_design(xt=xt1, groupsize=20, a=[[2,3,4], [4,5,6]])

        design_11 = create_design(xt=Matrix([0, 1, 2, 4, 6, 8, 24]), groupsize=50, 
                                    a=Matrix([70, 100], datanam=["WT", "DOSE"]))
        
        design_12 = create_design(xt=Matrix([0, 1, 2, 4, 6, 8, 24]), groupsize=50, 
                                    a=Matrix([70, 100], datanam=["WT", "DOSE"]), m=2)

        design_13 = create_design(xt=Matrix([0, 1, 2, 4, 6, 8, 24]), groupsize=50, 
                                    a=Matrix([[70, 100], [90, 200, 45]], datanam=[["WT", "DOSE"], ["DOSE", "WT", "AGE"]]), m=2)

        # design_14 same as design_13

        design_15 = create_design(xt=xt4, groupsize=np.array([50, 20]), 
                                    a=Matrix([[2, 3, 4], [4, 5, 6]], datanam=[["DOSE", "WT", "AGE"], [None, None, None]]))
if __name__ == '__main__':
    unittest.main()