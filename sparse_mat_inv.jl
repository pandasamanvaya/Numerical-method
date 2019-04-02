using LinearAlgebra

function gen_diag_mat(U)

	l = length(U[:,1])
	D = zeros(l,l)
	for i = 1:l
		D[i,i] = U[i,i]
		U[i,i] = 1
		for j = i+1:l
			U[i,j] = U[i,j]/D[i,i]
		end
	end

	return D, U 
end

function diag_inv(D)

	for i = 1:length(D[:,1])
		D[i,i] = 1/D[i,i]
	end

	return D
end

function sum_prod(A, B, i, j)
	sum = 0
	if i <= j
		for k = i+1:length(A[:,1])
			sum = sum + A[i,k]*B[k,j]
		end
	else
		for k = j+1:length(A[:,1])
			sum += A[i,k]*B[k,j]
		end
	end

	return sum
end

function sparse_mat_inv(A)
	F = lu(A)
	L = F.L
	D, U = gen_diag_mat(F.U)
	D_inv = diag_inv(D)
	Z = Matrix{Float64}(I, length(A[:,1]), length(A[:,1]))

	for i = length(A[:,1]):-1:1

		#Diagonal
		Z[i,i] = D_inv[i,i] - sum_prod(U,Z,i,i)

		#Upper triangle
		for j = length(A[1,:]):-1:i+1
			if A[i,j] != 0
				Z[i,j] = - sum_prod(U,Z,i,j)
			end
		end

		#Lower triangle
		for j = i-1:-1:1
			if A[i,j] != 0
				Z[i,j] = - sum_prod(Z,L,i,j)
			end
		end

	end
	return Z

end

#Test case
A = [6   96  102   124   -83;
  -9  -22   87    -6   103;
 -29   42   39   -34  -112;
   5   90   67   -60   118;
 113  -48   21  -104   -38
]

Z = sparse_mat_inv(A)
for i in 1:length(Z[:,1])
	println(Z[i,:])
end
