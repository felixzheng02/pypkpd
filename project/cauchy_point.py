"""


Author: Caiya Zhang, Yuchen Zheng
"""


import numpy as np
import math
from project.zeros import zeros


def cauchy_point(x, l, u, g, B):
    n = len(x)
    t = zeros(n, 1)
    d = zeros(n, 1)
    for i in range(0, n):
        if g[i] < 0:
            t[i] = (x[i] - u[i])/g[i]
        if g[i] > 0:
            t[i] = (x[i] - l[i])/g[i]
        if g[i] == 0:
            t[i] = math.inf
        if t[i] == 0:
            d[i] = 0
        else:
            d[i] = -g[i]
    ts_tmp, idx = np.unique(np.array(t), return_index=True) 
    ts = ts_tmp[idx.argsort()] 
    ts = ts[ts != 0]
    ts = np.sort(ts)

    df = np.cross(np.transpose(g), d)
    d2f = np.cross(np.cross(np.transpose(d), B), d)
    dt_min = -df/d2f
    t_old = 0
    i = 1
    z = zeros(n, 1)
    while (i <= len(ts) and dt_min >= (ts[i] - t_old)):
        ind = ts[i] < t
        d[~ind] = 0
        z = z + (ts[i] - t_old) * d
        df = np.cross(np.transpose(g), d) + np.cross(np.cross(np.transpose(d), B), z)
        d2f = np.cross(np.cross(np.transpose(d), B), d)
        dt_min = -df/d2f
        t_old = ts[i]
        i += 1
    dt_min = np.max(np.nanmax(dt_min), 0)
    t_old += dt_min
    x_cp = x - t_old * g
    temp = x - t * g
    x_cp[t_old > t] = temp[t_old > t]

    return x_cp