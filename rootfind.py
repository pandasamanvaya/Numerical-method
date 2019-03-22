#Biscetion Method#

def bisect (f, ainit, binit, err=1e-5, maxiter=50) :
    if f(ainit) * f(binit) > 0 :
        print 'Bracketing incorrect: change initial values'
        return 0
    
    if abs(f(ainit)) <= err :
        return (ainit)
    if abs(f(binit)) <= err :
        return (binit)
    
    i = 0
    xcur = binit
    ai = ainit
    bi = binit
    while i < maxiter :
        xcur = xcur + (ai - bi) / 2.0
        if abs(f(xcur)) <= err :      # Good solution found
            return xcur
        print 'x_', i, ': ', xcur, '  Error: ', f(xcur)
        if f(xcur) * f(ai) < 0 :   # Interval is [ai, xcur]
            bi = xcur
        elif f(xcur) * f(bi) < 0 : # Interval is [xcur, bi]
            ai = bi                # Reverse the limits: why?
            bi = xcur
        i = i + 1
    # End of while
    
    if i >= maxiter :
        print 'Good solution not found; Maximum number of iterations reached'
    
    return xcur

#End#

#Newton Method#
def newton(f, df, x, err = 1e-5, maxiter = 50):
    i = 0

    while i < maxiter:
        if abs(f(x)) == 0 :
            print 'Error = ',f(x)
            return x
        elif abs(f(x)) <= err:       #Solution Found
            print 'Error = ',f(x)
            return x

        print 'x_',i ,'= ',x ,'Error = ', f(x)
        
        x = x - f(x)/df(x)
        if abs(f(x)/df(x)) <= 1e-15:         #Divison by 0
            print 'Derivative is getting smaller than 1e-15'
            return x
        
        i = i + 1

        if i >= maxiter:
            print' Could not find solution in ', maxiter, 'iterations'
#End#

#Secant Method#
def secant(f, x1, x2, err = 1e-5, maxiter = 50):
    i = 0
    
    if f(x1) == 0:
        print i
        return x1
    elif abs(f(x1)) <= err:
        print i
        return x1
    
       
    while i < maxiter:
        if f(x2) == 0:
            return x2
        elif abs(f(x2)) <= err:      #Solution Found
            return x2

        if x2 == x1:
            print 'Both numbers are equal'
            if i != 0:
                return x2
            break

        der = (f(x2) - f(x1)) / (x2 - x1)
        if der == 0:            #Divison by 0
            print ' Derivative is very small'
            return x2
            
        print'x_',i, '=',x2, 'Error = ', f(x2)
        x3 = x2 - f(x2)/der
       
        x1 = x2
        x2 = x3

        i = i + 1

        if i == maxiter:
            print 'Could find solution in', maxiter,' iterations'
