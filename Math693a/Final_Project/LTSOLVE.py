def LTSOLVE(n,y,L):

	from numpy import zeros

	x = zeros(y.shape)

	x[n-1] = y[n-1]/L[n-1][n-1]

	for i in range(n-1,0,-1):

		x[i] = (y[i] - sum([L[j][i] * x[j] for j in range(i+1, n)]))/L[i][i]

	return x