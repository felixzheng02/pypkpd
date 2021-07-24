"""

## Function similar to that in matlab
## Author: Caiya Zhang
"""


import numpy as np
from project.feval import feval


def compute_step(stx,fx,dx,sty,fy,dy,stp,fp,dp,bracket,stpmin,stpmax):
    sgnd = dp*(dx/abs(dx))
    #First case: A higher function value. The minimum is bracketed. 
    #If the cubic step is closer to stx than the quadratic step, the 
    #cubic step is taken, otherwise the average of the cubic and 
    #quadratic steps is taken.
    if fp > fx:
        theta = 3*(fx-fp)/(stp-stx)+dx+dp
        s = max(np.array([abs(theta), abs(dx), abs(dp)]))
        gamma = s*np.sqrt((theta/s)^2-(dx/s)*(dp/s))
        if stp < stx:
            gamma=-gamma
        
        p = (gamma-dx) + theta
        q = ((gamma-dx) + gamma) + dp
        r = p/q
        stpc = stx + r*(stp-stx)
        stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp-stx)
        if abs(stpc-stx) < abs(stpq-stx):
            stpf = stpc
        else:
            stpf = stpc + (stpq-stpc)/2
        
        bracket = True
        #Second case: A lower function value and derivatives of opposite 
        #sign. The minimum is bracketed. If the cubic step is farther from
        #stp than the secant step, the cubic step is taken, otherwise the
        #secant step is taken.
    elif sgnd < 0:
        
        theta = 3*(fx - fp)/(stp - stx) + dx + dp
        s = max(np.array([abs(theta), abs(dx), abs(dp)]))
        gamma = s*np.sqrt((theta/s)^2 - (dx/s)*(dp/s))
        if stp > stx:
            gamma = -gamma
        
        p = (gamma - dp) + theta
        q = ((gamma - dp) + gamma) + dx
        r = p/q
        stpc = stp + r*(stx - stp)
        stpq = stp + (dp/(dp - dx))*(stx - stp)
        if abs(stpc-stp) > abs(stpq-stp):
            stpf = stpc
        else:
            stpf = stpq
        
        bracket = True
        #Third case: A lower function value, derivatives of the same sign,
        # and the magnitude of the derivative decreases.
    elif abs(dp) < abs(dx):
        #The cubic step is computed only if the cubic tends to infinity
        #in the direction of the step or if the minimum of the cubic
        #is beyond stp. Otherwise the cubic step is defined to be the
        #secant step.
        theta = 3*(fx - fp)/(stp - stx) + dx + dp
        s = max(np.array([abs(theta),abs(dx),abs(dp)]))
        
        #The case gamma = 0 only arises if the cubic does not tend
        #to infinity in the direction of the step.
        
        gamma = s*np.sqrt(max(0, (theta/s)^2-(dx/s)*(dp/s)))
        if stp > stx:
            gamma = -gamma
        
        p = (gamma - dp) + theta
        q = (gamma + (dx - dp)) + gamma
        r = p/q
        if r < 0 and gamma != 0:
            stpc = stp + r*(stx - stp)
        elif stp > stx:
            stpc = stpmax
        else:
            stpc = stpmin
        
        stpq = stp + (dp/(dp - dx))*(stx - stp)
        
        if bracket is True:
        
            #A minimizer has been bracketed. If the cubic step is
            #closer to stp than the secant step, the cubic step is
            #taken, otherwise the secant step is taken.
            
            if abs(stpc-stp) < abs(stpq-stp):
                stpf = stpc
            else:
                stpf = stpq
            
            if stp > stx:
                stpf = min(stp + 0.66*(sty-stp), stpf)
            else:
                stpf = max(stp + 0.66*(sty-stp), stpf)
            
        else:
            #A minimizer has not been bracketed. If the cubic step is
            #farther from stp than the secant step, the cubic step is
            #taken, otherwise the secant step is taken.
            
            if abs(stpc-stp) > abs(stpq-stp):
                stpf = stpc
            else:
                stpf = stpq
        
            stpf = min(stpmax, stpf)
            stpf = max(stpmin, stpf)
        
        #Fourth case: A lower function value, derivatives of the same sign, 
        #and the magnitude of the derivative does not decrease. If the
        #minimum is not bracketed, the step is either stpmin or stpmax,
        #otherwise the cubic step is taken.
    else:
        if bracket is True:
            theta = 3*(fp - fy)/(sty - stp) + dy + dp
            s = max(np.array([abs(theta), abs(dy), abs(dp)]))
            gamma = s*np.sqrt((theta/s)^2 - (dy/s)*(dp/s))
            if stp > sty:
                gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
        elif stp > stx:
            stpf = stpmax
        else:
            stpf = stpmin
        
    
    if fp > fx:
        sty = stp
        fy = fp
        dy = dp
    else:
        if sgnd < 0:
            sty = stx
            fy = fx
            dy = dx
        stx = stp
        fx = fp
        dx = dp
    
    stp = stpf
    return {"stx": stx, "fx":fx, "dx": dx, "sty": sty, "fy": fy, "dy": dy, "stp": stp, "bracket": bracket}


def line_search(f_name, f_options, l, u, x, f, g, d, options={}):
    #determine maximum step size
    
    if not any(options.keys()) is "ftol": 
        options["ftol"] = 1e-3
    
    if not any(options.keys()) is "gtol":
        options["gtol"] = 0.9
    
    if not any(options.keys()) is "xtol":
        options["xtol"] = 0.1
    
    
    xtol = options["xtol"]
    gtol = options["gtol"]
    ftol = options["ftol"]
    
    fixed = (x <= l or x >= u)
    stpmx = np.Infinity
    temp1 = np.Infinity
    for i in all(fixed == False):
        dk = d[i]
        if dk < 0:
            temp2 = l[i] - x[i]
            if temp2 >= 0:
                temp1 = 0
            else:
                if dk*stpmx < temp2: 
                    temp1=temp2/dk 
        else:
            temp2 = u[i] - x[i]
            if temp2 <= 0:
                temp1 = 0
            elif dk*stpmx is not np.nan:
                if dk*stpmx > temp2: 
                    temp1 = temp2/dk
    
        stpmx = min(temp1, stpmx)
    
    stp = min(1/np.linalg.norm(d, ord=np.Infinity), stpmx)
    stpmin = 0
    ##calc directional derivative
    gp = np.matmul(np.transpose(g), d)
    
    xtrapl = 1.1
    xtrapu = 4
    bracket = False
    stage = 1
    finit = f
    ginit = gp
    gtest = ftol*ginit
    width = stpmx - stpmin
    width1 = width/0.5
    stx = 0
    fx = finit
    gx = ginit
    sty = 0
    fy = finit
    gy = ginit
    stmin = 0
    stmax = stp + xtrapu*stp
    while True:
        f_options[0] = x + stp*d
        returnArgs = feval(f_name, f_options)
        f = returnArgs[0]
        g = returnArgs[1]
        gp = np.matmul(np.transpose(g), d)
        ftest = finit+stp*gtest
        if stage == 1 and f <= ftest and gp >= 0:
            stage = 2
        
        if f <= ftest and abs(gp) <= gtol*(-ginit):
            break
        
        if bracket is True and (stp <= stmin or stp >= stmax):
            break
        
        if bracket is True and stmax-stmin < xtol*stmax:
            break
        
        if stp == stpmx and f <= ftest and gp <= gtest:
            break
        
        if stp == stpmin and (f > ftest or gp >= gtest):
            break
        
        if stage == 1 and f <= fx and f > ftest:
            fm = f - stp*gtest
            fxm = fx - stx*gtest
            fym = fy - sty*gtest
            gm = gp - gtest
            gxm = gx - gtest
            gym = gy - gtest
            returnArgs = compute_step(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,bracket,stmin,stmax) 
            stx = returnArgs[0]
            fx = returnArgs[1]
            dx = returnArgs[2]
            sty = returnArgs[3]
            fy = returnArgs[4]
            dy = returnArgs[5]
            stp = returnArgs[6]
            bracket = returnArgs[7]
            fx = fxm + stx*gtest
            fy = fym + sty*gtest
            gx = gxm + gtest
            gy = gym + gtest
        else:
            returnArgs = compute_step(stx,fx,gx,sty,fy,gy,stp,f,gp,bracket,stmin,stmax) 
            stx = returnArgs[0]
            fx = returnArgs[1]
            dx = returnArgs[2]
            sty = returnArgs[3]
            fy = returnArgs[4]
            dy = returnArgs[5]
            stp = returnArgs[6]
            bracket = returnArgs[7]
        
        #Decide if a bisection step is needed.
        if bracket is True:
            if abs(sty-stx) > 0.66*width1:
                stp = stx + 0.5*(sty-stx)
            
            width1 = width
            width = abs(sty-stx)
        
        if bracket is True:
            stmin = min(stx, sty)
            stmax = max(stx, sty)
        else:
            stmin = stp + xtrapl*(stp-stx)
            stmax = stp + xtrapu*(stp-stx)
        
        stp = max(stp, stpmin)
        stp = min(stp, stpmx)
        if bracket and (stp <= stmin or stp >= stmax) or (bracket is True and (stmax-stmin) < xtol*stmax):
            stp = stx
        
        
    x = x + stp*d
    return {"f": f, "g": g, "x": x} 
    
