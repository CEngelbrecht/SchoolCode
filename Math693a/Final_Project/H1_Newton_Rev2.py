def H1_backtracking_and_fplus(n,x_c,f_c,g_c,s_N):

	from rosenbrock_2Nd_translated import rosenbrock
	import numpy as np

	#~~~~~~~~Backtracking~~~~~~~~~~~~~~~~~#

	alpha_bar = 1
	rho = 0.5
	c = 1E-4

	alpha = alpha_bar

	#function_value = eval_function(x_k + (p_k * alpha))

	function_value = rosenbrock(x_c + (s_N*alpha),0)

	test_value = rosenbrock((x_c),0) + (c * alpha * np.dot(np.transpose(s_N),g_c))


	while function_value > test_value: 

		alpha = rho * alpha

		function_value = rosenbrock(x_c + (s_N*alpha),0)

		test_value = rosenbrock((x_c),0) + (c * alpha * np.dot(np.transpose(s_N),g_c))

	#~~~~~~~~End Backtracking~~~~~~~~~~~~~~~~~#

	alpha_k = alpha #This is the alpha to use for this step in the process

	x_c = x_c + alpha_k*s_N #Iterating with step size gained from backtracking 

	retcode = 0 

	x_plus = x_c

	f_plus = rosenbrock(x_c,0)

	maxtaken = 0 

	return retcode,x_plus,f_plus,maxtaken