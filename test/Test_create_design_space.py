import numpy as np
from project.create_design import create_design
from project.create_design_space import create_design_space

design_1 = create_design(xt=[np.array([1,2,3,4,5]),
                             np.array([1,2,3,4])],
                         groupsize=np.array([50, 20]),
                         a=[np.array([70, 1000]),
                            np.array([35, 1000])]
                        #  a=[np.array([{"WT": 70, "DOSE": 1000}]),
                        #     np.array([{"DOSE": 1000, "WT": 35}])]
)

ds_1 = create_design_space(design_1)

ds_1_a = create_design_space(design_1, our_zero = 1e-5)

ds_2 = create_design_space(design_1, maxni=10, maxxt=10, minxt=0)

ds_3 = create_design_space(design_1,maxni=10, mingroupsize=20, maxxt=10, minxt=0)

ds_4 = create_design_space(design_1, maxa=np.array([100,2000]))

ds_5 = create_design_space(design_1, mina=np.array([10,20]))

design_2 = create_design(xt=[np.array([1,2,3,4,5]),
                             np.array([1,2,3,4])],
                         groupsize=np.array([50,20]),
                         a = [np.array([70, 1000]),
                              np.array([35, 1000])],
                         x = [np.array([1, 100]),
                              np.array([2, 200])]
)

ds_6 = create_design_space(design_2) 

ds_7 = create_design_space(design_2,
                           x_space=[np.array([1,2]),
                                    np.array(list(range(100,420,20)))])

ds_8 = create_design_space(design_2,
                           x_space=[np.array([1,2]),
                                    np.array(list(range(100,420,20)))],
                           grouped_xt=np.array([1,2,3,4,5])) 

ds_9 = create_design_space(design_2,
                           x_space=[np.array([1,2]),
                                    np.array(list(range(100,420,20)))],
                           use_grouped_xt=True) 

design_3 = create_design(xt=[np.array([1,2,3,4,5]),
                             np.array([1,2,3,4])],
                         groupsize=np.array([50,20]),
                         a = [np.array([35, 1000])],
                         x = [np.array([1, 100])]
)

ds_10 = create_design_space(design_3,
                           x_space=[np.array([1,2]),
                                    np.array(list(range(100,420,20)))],
                           use_grouped_a=True)

ds_11 = create_design_space(design_2,
                           x_space=[np.array([1,2]),
                                    np.array(list(range(100,420,20)))],
                           grouped_a=[np.array([1,2]), np.array([3,2])])

ds_12 = create_design_space(design_3,
                           x_space=[np.array([1,2]),
                                    np.array(list(range(100,420,20)))],
                           use_grouped_x=True)

ds_13 = create_design_space(design_3,
                           x_space=[np.array([1,2]),
                                    np.array(list(range(100,420,20)))],
                           grouped_x=[np.array([1,2]), np.array([3,2])])


# seq_1 = 1:10
# ds_14 = create_design_space(design_1,maxxt=10,minxt=0,
#                              xt_space = list(seq_1,seq_1,seq_1,seq_1,seq_1))
# ds_15 = create_design_space(design_1,maxxt=10,minxt=0,xt_space = list(seq_1))

# possible_values = as.matrix(cbind(list(0:10),list(0:10),list(0:10),list(0:20),list(0:20)))
# xt_space = as.matrix(rbind(possible_values,possible_values))

# ds_16 = create_design_space(design_1,maxxt=10,minxt=0,xt_space = xt_space)

# ds_17 = create_design_space(design_1,a_space = list(1:100,seq(1000,100000,by=1000)))

