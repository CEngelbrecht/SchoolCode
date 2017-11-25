def MODELHESS(n,machineps, H_c): 
	'''Need a function description'''

	import numpy as np
	from CHOLDECOMP import CHOLDECOMP

	sqrteps = np.sqrt(machineps)

	maxdiag = max([H_c[i][i] for i in range(len(H_c))])
	mindiag = min([H_c[i][i] for i in range(len(H_c))])
	maxposdiag = max(0,maxdiag)
	
	if mindiag <= sqrteps * maxposdiag: 

		mu = 2.0 * (maxposdiag - mindiag) * sqrteps - mindiag
		maxdiag = maxdiag + mu

	else: 

		mu = 0 

	#finding max off diagonal
	maxoff = 0

	for i in range(len(H_c)):
		for j in range(len(H_c)):
			if i == j:
				pass
			else: 
				if abs(H_c[i][j]) > maxoff: 
					maxoff = abs(H_c[i][j])

	if maxoff * (1.0 + 2.0*sqrteps) > maxdiag: 

		mu = mu + (maxoff - maxdiag) + 2.0 * sqrteps * maxoff
		maxdiag = maxoff * (1 + 2.0 * sqrteps)

	#9 
	if maxdiag == 0: #This means the Hessian is zero ?

		mu = 1.0
		maxdiag = 1.0 

	#10
	if mu > 0: 

		for i in range(len(H_c)):

			H_c[i][i] += mu

	#11 
	maxoffl = np.sqrt(max(maxdiag,(maxoff/n)))

	#12
	L,maxadd = CHOLDECOMP(n, H_c, maxoffl, machineps) #just does numpy cholesky factorization

	#13

	if maxadd > 0: 
		#do this later
		pass

	#13.8

	L, maxadd = CHOLDECOMP(n, H_c, 0, machineps)

	return L 





