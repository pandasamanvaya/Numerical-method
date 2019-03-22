import numpy as np
import eqsolve as eq

def simple(X, Y):

	sum_xx = (X*X).sum()
	sum_xy = (X*Y).sum()
	sum_x = X.sum()
	sum_y = Y.sum()

	N = len(X)

	A = np.array(([sum_xx, sum_x],
		[sum_x, N]))

	B = np.array(([sum_xy],
		[sum_y]))

	return(A, B)

def general(X,Y, n):
	sol = np.zeros((n,n))
	V = np.vander(X,n)
	
	l = len(X);
	for i in range(0,l):
		for j in range(0,n/2):
			k = V[i][j]
			V[i][j] = V[i][n-1-j]
			V[i][n-1-j] = k
	
	V_t = V.transpose()

	A = V_t.dot(V)
	B = V_t.dot(Y)
	
	sol = eq.gauss(A,B)
	return sol
