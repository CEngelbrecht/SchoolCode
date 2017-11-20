import numpy as np 
from rosenbrock_2Nd_translated import rosenbrock
import matplotlib.pyplot as plt 
import time 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def backtracking_search(p_k,x_k, grad): 
	#Note, the curvature condition for BFGS won't be satisfied here.

	inner_counter = 1.0

	alpha_bar = 1
	rho = 0.5
	c = 1E-4

	alpha = alpha_bar

	test_val1 = rosenbrock((x_k + (alpha * p_k)),0) #function evaluated at x_k + alpha * p_k 

	test_val2 = rosenbrock(x_k,0) + (c * alpha * np.transpose(p_k).dot(grad))

	while test_val1 > test_val2: 

		inner_counter += 1 
		
		alpha = rho * alpha

		test_val1 = rosenbrock((x_k + (alpha * p_k)),0)

		test_val2 = rosenbrock(x_k,0) + (c * alpha * np.transpose(p_k).dot(grad))

	return alpha,inner_counter

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def linesearch_H2(p_k,x_k,grad):
	pass


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def BFGS(): 

	'''Computes a pseudo Hessian'''

	alpha_list = [] 
	function_values = []
	grad_norm_list = []
	criteria_list = []

	inner_counter = 0

	time_start = time.time()

	x_0 = rosenbrock(0,-1) #first argument doesn't matter, -1 return the initial conditions 

	H_0 = np.identity(len(x_0))
		
	epsilon = 1E-8 #threshold 

	x_k = x_0 #initial conditions 
	H_k = H_0 

	grad = rosenbrock(x_k,1)
	grad_norm = np.linalg.norm(grad)

	I = np.identity(len(x_k)) #Initial guess to the inverse Hessian

	function_values.append(rosenbrock(x_k,0))
	grad_norm_list.append(grad_norm)

	inner_counter = 0.0
	outer_counter = 1.0

					#BFGS algorithm 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	while grad_norm > epsilon: 

		p_k = np.dot(-H_k,grad)

		alpha_k,temp_inner = backtracking_search(p_k,x_k,grad) #return step length

		x_kplusone = x_k + alpha_k * p_k #take step in p_k direction 

		s_k = x_kplusone - x_k

		y_k = rosenbrock(x_kplusone, 1) - grad

		rho_k = 1.0/np.dot(np.transpose(y_k),s_k)

		H_kplusone = (I - rho_k * np.dot(s_k,np.transpose(y_k))) * H_k * (I - rho_k * np.dot(y_k,np.transpose(s_k))) + rho_k * np.dot(s_k,np.transpose(s_k)) #BFGS update

		H_k = H_kplusone 

		x_k = x_kplusone 

		grad = rosenbrock(x_k,1) #update gradient and its norm 

		grad_norm = np.linalg.norm(grad)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


						#Convergence Criteria (might not work for backtracking)
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		criteria_list.append(np.linalg.norm((H_k - rosenbrock(x_k,2)) * p_k) / np.linalg.norm(p_k))
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


					#book keeping 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		grad_norm_list.append(grad_norm) #book keeping
		function_values.append(rosenbrock(x_k,0))
		outer_counter += 1
		inner_counter += temp_inner
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	#End BFGS
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	time_stop = time.time()
	function_values = [i[0][0] for i in function_values] #extract floats from a list of arrays 

	print("Done in {} seconds".format((time_stop - time_start)))

	plt.semilogy(function_values, label = "Objective value")
	plt.semilogy(grad_norm_list, label = "Norm of the gradient")
	plt.grid()
	plt.legend(loc = 'upper right')
	plt.xlabel("Iteration number")
	plt.title('BFGS Backtracking Linesearch Results \n Outer Iterations: {} Inner Iterations: {}'.format(outer_counter, inner_counter))
	plt.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Full_Newton():
	'''Explicitly computes the inverse of the full Hessian for p_k'''

	alpha_list = [] 
	function_values = []
	grad_norm_list = []

	time_start = time.time()

	outer_counter = 1 
	inner_counter = 0

	epsilon = 1E-8 #threshold 

	x_0 = rosenbrock(0,-1) #initial condition

	x_k = x_0

	function_values.append(rosenbrock(x_k,0))
 
	grad = rosenbrock(x_0,1) #gradient at starting point

	hessian = rosenbrock(x_0,2) #Hessian at starting point

	grad_norm = np.linalg.norm(grad) #norm of the gradient 

					#Full Newton algorithm 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	while grad_norm > epsilon: #Same condition as in BFGS

		p_k = (-1) * np.dot(np.linalg.inv(hessian),grad)

		alpha_k,temp_inner = backtracking_search(p_k,x_k,grad)

		x_kplusone = x_k + alpha_k * p_k

		grad = rosenbrock(x_kplusone,1)

		grad_norm = np.linalg.norm(grad)

		hessian = rosenbrock(x_kplusone,2)

		x_k = x_kplusone

		#				book keeping 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		outer_counter += 1 
		inner_counter += temp_inner
		function_values.append(rosenbrock(x_k,0))
		grad_norm_list.append(grad_norm)
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				#End Newton
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	time_stop = time.time()
	function_values = [i[0][0] for i in function_values] #extract floats from a list of arrays 

	print("Done in {} seconds".format((time_stop - time_start)))

	plt.semilogy(function_values, label = "Objective value")
	plt.semilogy(grad_norm_list, label = "Norm of the gradient")
	plt.grid()
	plt.legend(loc = 'upper right')
	plt.xlabel("Iteration number")
	plt.title('Full Newton Backtracking Linesearch Results \n Outer Iterations: {} Inner Iterations: {}'.format(outer_counter, inner_counter))
	plt.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

	BFGS()
	#Full_Newton()