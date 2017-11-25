def LSOLVE(n,b,L): 

	from numpy import zeros 

	# b = g_c

	y = zeros(b.shape)

	y[0] = b[0]/(L[0][0])

	for i in range(1,n):

		y[i] = (b[i] - sum([L[i][j] * y[j] for j in range(1,i-1)]))/L[i][i]

	return y 

