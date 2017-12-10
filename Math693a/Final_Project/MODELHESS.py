def MODELHESS(n,machineps, H_c): 
	'''Find mu necessary to make H + mu*I safely positive definite, and calculate Cholesky Decomposition of LLT of H'''

	from numpy import sqrt
	from CHOLDECOMP import CHOLDECOMP
	from numpy.linalg import cholesky 

	sqrteps = sqrt(machineps)

	maxdiag = max([H_c[i][i] for i in range(0,n)])
	mindiag = min([H_c[i][i] for i in range(0,n)])
	maxposdiag = max(0,maxdiag)
	
	if mindiag <= sqrteps * maxposdiag: 

		mu = 2.0 * (maxposdiag - mindiag) * sqrteps - mindiag
		maxdiag = maxdiag + mu

		print("mu = {}".format(mu))

	else: 

		mu = 0 

	#finding max off diagonal
	maxoff = 0

	for i in range(0,n):
		for j in range(0,n):
			if i == j: #don't look at the diagonal elements 
				pass
			else: 
				if abs(H_c[i][j]) > maxoff: 
					maxoff = abs(H_c[i][j]) #assign new value to maxoff

	if maxoff * (1.0 + 2.0*sqrteps) > maxdiag: 

		mu = mu + (maxoff - maxdiag) + 2.0 * sqrteps * maxoff
		maxdiag = maxoff * (1 + 2.0 * sqrteps)

	#9 
	if maxdiag == 0: #This means the Hessian is zero ?

		mu = 1.0
		maxdiag = 1.0 

	#10
	if mu > 0: 

		for i in range(0,n):

			H_c[i][i] += mu #add mu to every element 

	#11 
	maxoffl = sqrt(max(maxdiag,(maxoff/n)))

	#12
	#print("Calling CHOLDECOMP from MODELHESS")
	#L,maxadd = CHOLDECOMP(n, H_c, maxoffl, machineps) #does cholesky factorization 
	L = cholesky(H_c) #DOING NUMPY's CHOLESKY DOES EXACTLY THE SAME AS CHOLDECOMP 12/6/2017
	maxadd = 0 
	#13

	if maxadd > 0: 
		#do this later
		print("Maxadd > 0")
		L, maxadd = CHOLDECOMP(n, H_c, 0, machineps)
		pass

	#13.8

	return L 