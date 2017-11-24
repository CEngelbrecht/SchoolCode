import numpy as np 

def UMSTOP0(n, x_0, f_0,g_0,gradtol): 

	if max(np.linalg.norm(g_0), np.linalg.norm(x_0)/np.linalg.norm(f_0))< 1E-3 * gradtol:  

		termcode = 1 
	else: 
		termcode = 0 

	return termcode 