from Powell_function import Powell_function
import numpy as np
import matplotlib.pyplot as plt 

from Rosenbrock_Driver import FN_rosenbrock as FN
from Rosenbrock_Driver import GRAD_rosenbrock as GRAD
from Rosenbrock_Driver import HESS_rosenbrock as HESS

from FDGRAD import FDGRAD
from MACHINEPS import MACHINEPS

eta = MACHINEPS()

function = "Rosenbrock"

if function == "Powell":
	
	x_0 = (3,-1,0,1)

	x_0 = np.reshape(np.array(x_0),(len(x_0),1)) #numpy array of (x_0)T

	x_c = x_0

	print x_c

	f_list = []

	f_value = Powell_function(x_c,0)

	g_c = Powell_function(x_c,1)

	hess = Powell_function(x_c,2)

	f_list.append(f_value)

	while np.linalg.norm(g_c) > 1E-8: 

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




elif function == "Rosenbrock": 

	n = 2 

	x_0 = (-1.2,1.0)

	x_0 = np.reshape(np.array(x_0),(len(x_0),1)) #numpy array of (x_0)T

	f_list = []

	f_value = FN(n,x_0)

	#g_c = GRAD(n,x_0)
	g_c = FDGRAD(n,x_0,f_value,FN,eta)

	hess = HESS(n,x_0)
	f_list.append(f_value)

	x_c = x_0


	while np.linalg.norm(g_c) > 1E-8: 

		s_N = -1 * np.dot(np.linalg.inv(hess),g_c)

	#~~~~~~~~Backtracking~~~~~~~~~~~~~~~~~#

		alpha_bar = 1
		rho = 0.5
		c = 1E-4

		alpha = alpha_bar

		function_value = FN(n,x_c + (s_N*alpha))

		test_value = FN(n,(x_c)) + (c * alpha * np.dot(np.transpose(s_N),g_c))

		while function_value > test_value: 

			alpha = rho * alpha

			function_value = FN(n,(x_c + (s_N*alpha)))

			test_value = FN(n,(x_c)) + (c * alpha * np.dot(np.transpose(s_N),g_c))

	#~~~~~~~~End Backtracking~~~~~~~~~~~~~~~~~#

		alpha_k = alpha #This is the alpha to use for this step in the process

		x_c = x_c + alpha_k*s_N #Iterating with step size gained from backtracking 

		f_value = FN(n,x_c)

		#g_c = GRAD(n,x_c)

		g_c = FDGRAD(n,x_c,f_value,FN,eta)

		hess = HESS(n,x_c)

		f_list.append(f_value)

	print x_c
	plt.semilogy([float(i) for i in f_list])
	plt.show()