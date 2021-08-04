
import time
from project.tictoc import tic
from project.tictoc import toc


tic()
toc()
#> Elapsed time is 0.0 seconds.
tic()
time.sleep(0.5)
toc()
#> Elapsed time is 0.5 seconds.
