import numpy as np
import time

#Gauss Elimination Method#
def gauss(A, B):
        #st_time = time.time()
        N = len(B)
        A_aug = gauss_augmat(A,B,N)
        A_aug = gauss_eliminate(A_aug, N)
        A_sol = gauss_substitute(A_aug, N)
        #time_calc(st_time)
        print 'Solutions are returned'
        return A_sol

def gauss_augmat(A , B, N):

    A_aug = np.zeros((N, N+1))
    for i in range(0, N):
        for j in range(0, N+1):
            if j <= N-1:
                A_aug[i][j] = A[i][j]
            else:
                A_aug[i][j] = B[i]

    return A_aug
def gauss_eliminate(A_aug, N):

    for i in range(0 , N-1):
        for j in range(i+1, N):
                fact = A_aug[j][i]/A_aug[i][i]

                A_aug[j] = A_aug[j] + ( -fact * A_aug[i])
                A_aug[j][i] = 0

    return A_aug

def gauss_substitute(A_aug, N):

    A_sol = np.zeros((N,1))
    for i in range(N-1, -1, -1):
        diff = A_aug[i][N]
        for j in range(N-1, i, -1):

            diff = diff - A_aug[i][j]*A_sol[j]

        A_sol[i] = diff/A_aug[i][i]

    return A_sol

#end of Gauss Elimination

#Jacobi Method#

def jacobi(A,B,x1,err = 1e-10):

    N = len(A)
    B = B.reshape(N,1)
    x1 = x1.reshape(N,1)

    D =  jacobi_diag(A, N)
    A_off = jacobi_off(A,D)
    D_inv = jacobi_inv(D,N)
    A_sol = jacobi_root(D_inv, B , A_off, x1, N, err)

    print 'Solutions are returned'
    return A_sol

def jacobi_diag(A, N):

    D = np.zeros((N,N))
    
    for i in range(0,N):
        for j in range(0,N):

            if i != j:
                D[i][j] = 0
            else:
                D[i][j] = A[i][j]
    return D

def jacobi_off(A, D):
    
    A_off = A - D

    return A_off

def jacobi_inv(D, N):

    D_inv = np.zeros((N,N))
    
    for i in range(0, N):
        for j in range(0, N):

            if i == j:
                D_inv[i][j] = 1/D[i][j]

            else:
                D_inv[i][j] = 0
    return D_inv
def jacobi_root(D_inv, B , A_off, x1, N, err):

    diff = np.zeros(N)
       
    for j in range(0, 100):
        x2 = D_inv.dot( B - A_off.dot(x1))

        diff = x2 - x1
        
        for i in range(0, N):
            val = False
            if abs(diff[i]) <= err:
                val = True
            else:
                break
        if val == True:
            print 'No.of iterations =', j+1
            return x2
        else:
            x1 = x2
        
          
    print 'Could not find after no.of iterations =', j
    return x2
#end of Jacobi

#L-U Decompostion#

def lu_decompose(A,B):

    N = len(A)
    (A_up,A_low) = lu_eliminate(A,N)
    A_y = lu_sub_low(A_low,B,N)
    A_sol = lu_sub_up(A_up,A_y,N)
    
    print 'Solutions are returned'
    return A_sol

def lu_eliminate(A,N):
    k = 0
    A_low = np.zeros((N,N))
    A_up = np.zeros((N,N))
    
    for i in range(0 , N):
         A_up[i] = A[i]

    for i in range(0 , N):
        for j in range(i+1, N):
            
            fact = A_up[j][i]/A_up[i][i]
            A_low[j][k] = fact
            A_up[j] = A_up[j] + ( -fact * A_up[i])
        
        k = k +1
        A_low[i][i] = 1
                
    return (A_up, A_low)

def lu_sub_low(A_low, B, N):

    A_y = np.zeros((N,1))

    for i in range(0, N):
        diff = B[i]
        for j in range(0, i):
        
            diff = diff - A_low[i][j]*A_y[j]
           
        A_y[i] = diff/A_low[i][i]

    return A_y

def lu_sub_up(A_up, A_y, N):

    A_sol = np.zeros((N,1))

    for i in range(N-1, -1, -1):
        diff = A_y[i]
        for j in range(N-1, i, -1):
        
            diff = diff - A_up[i][j]*A_sol[j]
           
        A_sol[i] = diff/A_up[i][i]

    return A_sol

#end of L-U Decomposition
 
#Gauss-Seidel# 

def gauss_seidel(A, B, x, err = 1e-10, maxiter = 50):
    #st_time = time.time()
    N = len(A)
    for k in range(0, maxiter):
        val = False

        sum1 = 0
        sum2 = 0
        
        for i in range(0, N):
            sum1 = 0
            sum2 = 0
            
            for j in range(0 , i):
                sum1 = sum1 + A[i][j] * x[j]
                
            for j in range(i+1, N):
                sum2 = sum2 + A[i][j] * x[j]
            
            sol = (B[i] - sum1 - sum2)/A[i][i]
            if abs(sol - x[i]) <= err:
                val = True
            else :
                x[i] = sol
		val = False
	
        if val == True:
            print 'No.of iterations =', k
            x = x.reshape(N,1)
            print 'Solutions are returned'
            #time_calc(st_time)
            return x
    print 'Found no solution in',k,'iterations'

#end of Gauss-Seidel

def time_calc(st_time):
        print 'Time taken by method = ',time.time()- st_time,'s'
        return
