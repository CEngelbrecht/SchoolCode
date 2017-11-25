def CHOLDECOMP(n, H_c, maxoffl, machineps):
	'''Currently just uses numpy'''

	from numpy.linalg import cholesky
	from numpy import sqrt

	minl = (machineps)**(1.0/4.0) * maxoffl

	if maxoffl == 0: 

		maxoffl = sqrt(max([H_c[i][i] for i in range(len(H_c))]))

	minl2 = sqrt(machineps) * maxoffl

	maxadd = 0 

	L = cholesky(H_c)

	return L,maxadd
	