def UMINCK (n,machineps,x_0,**kwargs):
	'''Checks values of algorithmic options and tolerances'''

	if 'fdigits' in kwargs:
		fdigits = kwargs['fdigits']

		if fdigits == -1:
			eta = machineps
		else: 
			eta = max(machineps,1*10**(-fdigits)) 
	
	if n < 1: 
		termcode = -1
		return termcode

	else:
		termcode = 0 


	return termcode,eta


