def HESS(n,x_k): 

	from rosenbrock_2Nd_translated import rosenbrock

	HESS = rosenbrock(x_k,2)

	return HESS 

