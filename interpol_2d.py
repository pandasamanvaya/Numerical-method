import numpy as np
import interpolation as int

def interpol(X,Y,Z,pt):

	n = len(X)
	A = np.zeros(n)
	B = np.zeros(n)
	L = np.zeros(n)
	C = np.zeros(n)
	i = 0
	j = 0
	k = 0
	while i < n and j < n:
		if X[j] == X[i]:
			A[j-i] = Y[j]
			B[j-i] = Z[j]
			j = j + 1
			flag = True
			if j == n:
				flag = False
		else:	
			flag = False
		
		if flag == False:
			print j,pt[:,1],A[:4],B[:4]
			L[k] = int.vander(pt[:,1],A[:4],B[:4])
			C[k] = X[i]
			i = j
			k = k + 1
			
	print pt[:,0],C[:k],L[:k]
	val = int.vander(pt[:,0],C[:k],L[:k])
	print val
	return 
