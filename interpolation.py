import numpy as np
import eqsolve as eq

def linear( pt, X, Y):

	n1 = len(X)
	n2 = len(pt)
	y = np.zeros((n2))
	for j in range (0, n2):
		if pt[j] <= X[0] or pt[j] > X[n1-1]:
			print 'One of values is not between',
			print X[0],' - ',X[n1-1]

			return

		for i in range (0,n1-1):
			if pt[j] >= X[i] and pt[j] < X[i+1]:
				break
	 
	  	ai = (Y[i+1] - Y[i]) / (X[i+1] - X[i])
	  	bi = Y[i]

	  	y[j] = (ai*(pt[j]-X[i])) + bi
	
	return y


def vander(pt, X, Y):

	y = 0
	n1 = len(X)
	n2 = len(pt)
	y = np.zeros((n2))

	A = np.zeros(n1)
	V = np.vander(X, n1)
                
        for i in range(0,n1):
                for j in range(0,n1/2):
                        k = V[i][j]
                        V[i][j] = V[i][n1-1-j]
                        V[i][n1-1-j] = k
        
        A = eq.gauss(V,Y)
        
        for j in range(0,n2):
                if pt[j] < X[0] or pt[j] > X[n1 - 1]:
                        print 'One of values is not between',
                        print X[0],' - ',X[n1-1]
		
                for i in range(0, n1):
			y[j] = (A[i]*(pt[j]**i)) + y[j]

	return y

def lagrange(pt , X, Y):

	n1 = len(X)
	n2 = len(pt)
	y = np.zeros(n2)
	
	for k in range(0,n2):
		if pt[k] < X[0] or pt[k] > X[n1 - 1]:
			print 'One of values is not between',
			print X[0],' - ',X[n1-1]
			return

		for i in range(0,n1):
			prod = 1.0
			for j in range(0,n1):

				if i != j :
					prod = (pt[k]-X[j])/(X[i]-X[j]) * prod
                        
                        y[k] = Y[i]*prod + y[k]
                        
	return y
