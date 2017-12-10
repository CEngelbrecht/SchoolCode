def CHOLSOLVE(n,g_c,L_c):
	'''Solve (L * LT)s = -g for s'''

	from LSOLVE import LSOLVE
	from LTSOLVE import LTSOLVE 

	s = LSOLVE(n,g_c,L_c)

	s = LTSOLVE(n,s,L_c)

	s = -1.0*s

	#print("s as returned by CHOLSOLVE: {}".format(s))
	return s 
	

