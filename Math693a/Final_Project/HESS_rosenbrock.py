def HESS_rosenbrock(n,x_k): 

	from rosenbrock_2Nd_translated import rosenbrock

	H_k = rosenbrock(x_k,2)

	return H_k

