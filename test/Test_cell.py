import numpy as np
from project.cell import cell


print(cell(3))
#>      [,1]      [,2]      [,3]     
#> [1,] numeric,0 numeric,0 numeric,0
#> [2,] numeric,0 numeric,0 numeric,0
#> [3,] numeric,0 numeric,0 numeric,0
print(cell(2,3))
#>      [,1]      [,2]      [,3]     
#> [1,] numeric,0 numeric,0 numeric,0
#> [2,] numeric,0 numeric,0 numeric,0
## define possible values of 2 categorical design variable
x_space = cell(1,2)
print(x_space)
"""
x_space[0,0] = np.arange(10,100+10,10).tolist()
x_space[1,2] = np.arange(10,300+10,10).tolist()
print(x_space)
"""
#>      [,1]       [,2]      
#> [1,] numeric,10 numeric,30x.space[1,1]
#> [[1]]
#>  [1]  10  20  30  40  50  60  70  80  90 100
#> x.space[1,2]
#> [[1]]
#>  [1]  10  20  30  40  50  60  70  80  90 100 110 120 130 140 150 160 170 180 190
#> [20] 200 210 220 230 240 250 260 270 280 290 300
#> 