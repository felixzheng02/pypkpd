import numpy as np
import project.create_design
import project.ones
from matpy.matrix import matrix


xt1 = [np.array([1, 2, 3]), np.array([1, 2, 3, 4])]


design_1 = project.create_design.create_design(xt=xt1, groupsize=20)
design_6 = project.create_design.create_design(xt=xt1, groupsize=20, model_switch=project.ones.ones(2, 4))
design_7 = project.create_design.create_design(xt=xt1, groupsize=20, a=matrix(np.array([2, 3, 4])))
design_8 = project.create_design.create_design(xt=xt1, groupsize=20, a=matrix(np.array([[2,3,4], [4,5,6]])))
design_9 = project.create_design.create_design(xt=xt1, groupsize=20, a=[matrix(np.array([2,3,4,6])), matrix(np.array([4,5,6]))])
design_10 = project.create_design.create_design(xt=xt1, groupsize=20, a=[matrix(np.array([2,3,4])), matrix(np.array([4,5,6]))])