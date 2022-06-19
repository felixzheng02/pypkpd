import sys
sys.path.append("/Users/felix/Documents/college/pypkpd")
import unittest
import numpy as np
import project.create_design
import project.ones
from matpy.matrix import Matrix


class TestCreateDesign(unittest.TestCase):

    def test_create_design_init(self):
        
        xt1 = [np.array([1, 2, 3]), np.array([1, 2, 3, 4])]
        xt4 = [np.array([1,2,3,4,5]), np.array([1,2,3,4])]
        xt2 = Matrix(np.array([[1,2,3,4], [1,2,3,4]]))
        xt3 = Matrix(np.array([1,2,3,4]))

        design_1 = project.create_design.create_design(xt=xt1, groupsize=20)
        design_2 = project.create_design.create_design(xt=xt4, groupsize=20)
        design_3 = project.create_design.create_design(xt=xt2, groupsize=20)

        #design_4 duplicates
        design_4 = project.create_design.create_design(xt=xt3, groupsize=20)

        design_5 = project.create_design.create_design(xt=xt3, groupsize=20,m=3)
        design_6 = project.create_design.create_design(xt=xt1, groupsize=20, model_switch=project.ones.ones(2, 4))
        design_7 = project.create_design.create_design(xt=xt1, groupsize=20, a=Matrix(np.array([2, 3, 4])))
        design_8 = project.create_design.create_design(xt=xt1, groupsize=20, a=Matrix(np.array([[2,3,4], [4,5,6]])))

        # design_9 warningas and "a" duplicates
        design_9 = project.create_design.create_design(xt=xt1, groupsize=20, a=[Matrix(np.array([2,3,4,6])), Matrix(np.array([4,5,6]))])

        design_10 = project.create_design.create_design(xt=xt1, groupsize=20, a=[Matrix(np.array([2,3,4])), Matrix(np.array([4,5,6]))])

        for key,values in design_9.items():
            if type(values) is Matrix:
                print(key, values.get_all_data())
            else:
                print(key, values)


    

if __name__ == '__main__':
    unittest.main()
    

