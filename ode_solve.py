import numpy as np
import eqsolve as eq

def ode_euler(f, tn, tdelta, x0 = 0, t0 = 0):

        i = t0
        xn = x0

        while (tn-i) > 1e-5:
                xn = xn + tdelta*f(i, xn)
                i = i + tdelta

        return i, xn

def ode_rk1(f, tn, tdelta, x0=0, t0=0):

        i = t0
        xn = x0

        while (tn-i) > 1e-5:
                k1 = tdelta * f(i, xn)
                k2 = tdelta * f(i + tdelta, xn + k1)

                xn = xn + 0.5*(k1 + k2)
                i = i + tdelta

        return i , xn

def ode_rk2_heun(f, tn, tdelta, x0 = 0, t0 = 0):

        i = t0
        xn = x0

        while (tn-i) > 1e-5:
                k1 = 0.5 * tdelta * f(i, xn)
                k2 = tdelta * f(i + 0.5*tdelta, xn + k1)

                xn = xn + k2
                i = i + tdelta

        return i, xn

def ode_rk2_midpoint(f, tn, tdelta, x0 = 0, t0 = 0):

        i = t0
        xn = x0

        while (tn-i) > 1e-5:
                k1 = tdelta * f(i, xn)
                k2 = tdelta * f(i + tdelta, xn + k1)

                xn = xn + 0.5*(k1 + k2)
                i = i + tdelta

        return i , xn

def ode_rk2_ralston(f, tn, tdelta, x0 = 0, t0 = 0):

        i = t0
        xn = x0

        while (tn-i) > 1e-5:
                k1 = tdelta * f(i, xn)
                k2 = tdelta * f(i + 0.75*tdelta, xn + 0.75*k1)

                xn = xn + 0.33333*(k1 + 2*k2)
                i = i + tdelta

        return i , xn

def ode_solve(f, tn, tdelta, x0 = 0, t0 = 0):

        i = t0
        xn = x0

        while (tn - i) > 1e-5:
                k1 = tdelta * f(i, xn)
                k2 = tdelta * f(i + 0.5*tdelta, xn + 0.5*k1)
                k3 = tdelta * f(i + 0.5*tdelta, xn + 0.5*k2)
                k4 = tdelta * f(i + tdelta, xn + k3)

                xn = xn + (1/6.0)*(k1 + 2*k2 + 2*k3 + k4)
                i = i + tdelta

        return i , xn
def ode_rk4_38(f, tn, tdelta, x0 = 0, t0 = 0):

        i = t0
        xn = x0

        while (tn-i) > 1e-5:
                k1 = tdelta * f(i, xn)
                k2 = tdelta * f(i + (1/3)*tdelta, xn + (1/3)*k1)
                k3 = tdelta * f(i + (2/3)*tdelta, xn + k2 - k1/3)
                k4 = tdelta * f(i + tdelta, xn + k1 - k2 + k3)

                xn = xn + (1/8.0)*(k1 + 3*k2 + 3*k3 + k4)
                i = i + tdelta

        return i , xn
def ode_fin_diff(f, k1, k2, y_k1, y_k2, n):

        A = np.zeros((n,n))
        B = np.zeros(n)
        y = np.zeros(n)

        h = (k2 - k1)/(n+1)

        B[0]= y_k1
        B[n-1] = y_k2

        for i in range (0,n):
                for j in range (0,n):
                        if abs(i-j) == 1 :
                                A[i][j] = -1
                        elif i==j:
                                A[i][j] = 2
        
        for i in range (0,n):
            B[i] = h**2 * f((i+1)*h) + B[i]
        
        
        y = eq.gauss(A,B)
        
        print"x_ 0  = ",k1 ," y_ 0 = ",y_k1
        for i in range(0,n):
                print "x_",i+1," = ",(i+1)*h, " y_",i+1," = ",y[i]
        print "x_",n+1," = ",k2 ," y_",n+1," = ",y_k2
        return 

def ode_shoot(f, b1, b2, h, x0, xn, y0, yn, err = 1e-5, maxiter = 50):

        F1 = calcF(f,x0,xn,y0,yn,b1,h)
        if abs(F1) <= err:
                b = b1

        F2 = calcF(f,x0,xn,y0,yn,b2,h)
        if abs(F2) <= err:
                b = b2

        for i in range(0,maxiter):
                der = (F2 - F1)/(b2 - b1)
                if der == 0:
                        print "Derivative is small"
                        break
                else:
                        b3 = b2 - F2/der
                F3 = calcF(f,x0,xn,y0,yn,b3,h)
                print "F = ",F3
                if abs(F3) <= err:
                        b = b3
                        break
                else:
                        F1 = F2
                        F2 = F3
                        b1 = b2
                        b2 = b3
        print "b = ",b
        print "i = ",i
        if(i == (maxiter-1)):
                print "Couldn't find solutions in",maxiter,"iterations."
        else:
                z0 = b
                j = x0
                print " Solutions are:-"
                print "x = ",j,"y = ",y0
                while (xn-j) > 1e-5:
                        y1 = y0 + h*f(x0, y0, z0)
                        z1 = z0*(h+1)
                        z0 = z1
                        y0 = y1

                        j = j + h
                        print "x = ",j,"y = ",y1
        return 


def calcF(f,x0,xn,y0,yn,b,h):
        z0 = b
        i = x0
        while (xn-i) > 1e-5:
                y1 = y0 + h*f(x0, y0, z0)
                z1 = z0*(h+1)
                z0 = z1
                y0 = y1
                i = i + h
        F = y1 - yn
        return F

