from Powell_function import Powell_function
import numpy as np
import matplotlib.pyplot as plt 

x_0 = (3,-1,0,1)

x_0 = np.reshape(np.array(x_0),(len(x_0),1)) #numpy array of (x_0)T

x_c = x_0

f_list = []

f_value = Powell_function(x_c,0)

g_c = Powell_function(x_c,1)

hess = Powell_function(x_c,2)

f_list.append(f_value)

while f_value > 1E-12: 

	s_N = -1 * np.dot(np.linalg.inv(hess),g_c)

#~~~~~~~~Backtracking~~~~~~~~~~~~~~~~~#

	alpha_bar = 1
	rho = 0.5
	c = 1E-4

	alpha = alpha_bar

	function_value = Powell_function(x_c + (s_N*alpha),0)

	test_value = Powell_function((x_c),0) + (c * alpha * np.dot(np.transpose(s_N),g_c))

	while function_value > test_value: 

		alpha = rho * alpha

		function_value = Powell_function(x_c + (s_N*alpha),0)

		test_value = Powell_function((x_c),0) + (c * alpha * np.dot(np.transpose(s_N),g_c))

#~~~~~~~~End Backtracking~~~~~~~~~~~~~~~~~#

	alpha_k = alpha #This is the alpha to use for this step in the process

	x_c = x_c + alpha_k*s_N #Iterating with step size gained from backtracking 

	f_value = Powell_function(x_c,0)

	g_c = Powell_function(x_c,1)

	hess = Powell_function(x_c,2)

	f_list.append(f_value)

print x_c
plt.semilogy([float(i) for i in f_list])
plt.show()
