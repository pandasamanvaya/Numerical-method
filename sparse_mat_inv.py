import scipy.linalg as linalg
import numpy as np

def recur_inv(Z):
	P, L, U = linalg.lu(Z)

	D = np.diag(np.diag(U))
	U /= np.diag(U)[:,None]

	D_inv = np.linalg.inv(D)
	U_inv = np.linalg.inv(U)

	err = 1
	delta = 100
	maxiter = 100
	while maxiter:
		Z = U_inv.dot(D_inv) + Z.dot(np.identity(len(Z)) - P.dot(L))
		maxiter -= 1

	return Z 

if __name__ == "__main__":
	
	n = int(input("Enter dimension : "))
	Z = [[int(i) for i in input().split()]for j in range(n)]
	Z = np.asarray(Z)
	print(Z)
	print(recur_inv(Z))