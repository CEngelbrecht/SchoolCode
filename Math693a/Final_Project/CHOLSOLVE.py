def CHOLSOLVE(n,g_c,L_c):

	from LSOLVE import LSOLVE
	from LTSOLVE import LTSOLVE 

	s = LSOLVE(n,g_c,L_c)

	s = LTSOLVE(n,s,L_c)

	s = -1*s

	return s 

