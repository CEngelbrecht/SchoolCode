def FN_Powell(n,x_c): 

	from Powell_function import Powell_function

	x_c = Powell_function(x_c,0)

	return x_c

def GRAD_Powell(n,x_c):

	from Powell_function import Powell_function

	g_c = Powell_function(x_c,1)

	return g_c

def HESS_Powell(n,x_c):

	from Powell_function import Powell_function

	H_c = Powell_function(x_c,2)

	return H_c