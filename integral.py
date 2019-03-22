from math import *

def int_all(func, a, b, h):
        mid = midpnt(func, a, b)
        trap = trapezoidal(func, a, b)
        sim = simpson(func, a, b)
        cmp_trap = comp_trapezoidal(func, a, b, h)
        cmp_sim = comp_simpson(func, a, b, h)

        print "Midpoint = ",mid
        print "Trapezoidal = ",trap
        print "Simpson  = ",sim
        print "Composite Trapezoidal = ",cmp_trap
        print "Composite Simpson = ", cmp_sim
        
        return

def midpnt(func, a, b):
        return((b - a) * func((a + b)/2.0))


def trapezoidal(func, a, b):
        return((b-a)/2.0 * (func(a) + func(b)))

def simpson(func, a, b):
        return((b-a)/6.0) * (func(a) + 4*func((a+b)/2.0) + func(b))

def comp_trapezoidal(func, a, b, h):

        x_i = a
        sum = 0.0
        while x_i < (b-h):

                x_i = x_i + h
                sum = 2.0 * func(x_i) + sum
        
        return (h/2) * (func(a) + sum + func(b))

def comp_simpson(func, a, b, h):
        x_i = a
        sum = 0.0
        while x_i < (b-h):

                x_i = x_i + h

                if trunc((x_i - a)/h) % 2 == 0:
                        sum = 2 * func(x_i) + sum
                else:
                        sum = 4 * func(x_i) + sum

                
        return (h/3.0) * (func(a) + sum + func(b))

