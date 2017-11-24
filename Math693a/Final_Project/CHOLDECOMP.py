def CHOLDECOMP(n, H_c, maxoffl, machineps):
	'''Currently just uses numpy'''

	from numpy.linalg import cholesky

	minl = (machineps)**(1.0/4.0)

	if maxoffl = 0: 

		maxoffl = np.sqrt(max([H_c[i][i] for i in range(len(H_c))]))

	minl2 = np.sqrt(machineps) * maxoffl

	maxadd = 0 

	L = cholesky(H_c)

	return L
	