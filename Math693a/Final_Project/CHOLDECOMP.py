def CHOLDECOMP(n, H_c, maxoffl, machineps):
	'''Solves for the cholesky decomposition of a matrix.
		Changelog: 25/11/2017 4:55 PM fixed indices on lines 24,38, now returns same for a 2x2 positive definite matrix as numpy.linalg.cholesky'''

	from numpy.linalg import cholesky
	from numpy import sqrt,zeros

	minl = (machineps)**(1.0/4.0) * maxoffl

	if maxoffl == 0: 

		maxoffl = sqrt(max([H_c[i][i] for i in range(len(H_c))]))

	minl2 = sqrt(machineps) * maxoffl

	maxadd = 0 

	#print("H_c = {} from CHOLDECOMP".format(H_c))

	L = zeros(H_c.shape)

	for j in range(0,n):

		L[j][j] = H_c[j][j] - sum([(L[j][i])**2 for i in range(0,j)])
		minljj = 0 

		for i in range((j+1),n):
			L[i][j] = H_c[j][i] - sum([L[i][k] * L[j][k] for k in range(0,j)]) #<- original 
			minljj = max(abs(L[i][j]),minljj)

		minljj = max((minljj/maxoffl),minl)

		if L[j][j] > minljj**2: 
			L[j][j] = sqrt(L[j][j]) #normal Cholesky iteration

		else: #Augment H[j][j]
			print("Augmenting H_c")

			if minljj < minl2: 
				minljj = minl2

			maxadd = max(maxadd,(minljj**2 - L[j][j]))
			L[j][j] = minljj
		for i in range(j+1,n):
			L[i][j] = L[i][j]/L[j][j]
	return L,maxadd
	