def FN_rosenbrock(n,x_c): 

	from rosenbrock_2Nd_translated import rosenbrock

	x_c = rosenbrock(x_c,0)

	return x_c

def GRAD_rosenbrock(n,x_c): 

	from rosenbrock_2Nd_translated import rosenbrock

	g_c = rosenbrock(x_c,1)

	return g_c


def HESS_rosenbrock(n,x_c): 

	from rosenbrock_2Nd_translated import rosenbrock

	H_c = rosenbrock(x_c,2)

	return H_c

