'''
Christian Engelbrecht
Fall 2017 Math 693a 
Homework 5
'''

import numpy as np 
from rosenbrock_2Nd_translated import rosenbrock
import matplotlib.pyplot as plt 
import time 

#						Start Backtracking Linesearch Function Definition
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def backtracking_search(p_k,x_k, grad): 
	'''Simple backtracking line search, a la Homework 1'''

	inner_counter = 1.0 #book keeping for tracking how many inner iterations there are 
	alpha_bar = 1 #standard values, always test 1 first 
	rho = 0.5
	c = 1E-4
	alpha = alpha_bar 

	test_val1 = rosenbrock((x_k + (alpha * p_k)),0) #function evaluated at x_k + alpha * p_k 
	test_val2 = rosenbrock(x_k,0) + (c * alpha * np.transpose(p_k).dot(grad)) #test condition 

	#					Start finding alpha
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	while test_val1 > test_val2: 

		inner_counter += 1 
		
		alpha = rho * alpha

		test_val1 = rosenbrock((x_k + (alpha * p_k)),0)

		test_val2 = rosenbrock(x_k,0) + (c * alpha * np.transpose(p_k).dot(grad))

	#					End finding alpha 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	return alpha,inner_counter

#						End Backtracking Linesearch Function Definition
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#						Start Definition for Strong Wolfe Search Function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def LS_Strong_Wolfe(p_k,x_k,grad):
	'''Homework 2 styled line search, with a simple zoom function'''	

	inner_counter = 1 

	phi = lambda alpha: rosenbrock((x_k + (alpha * p_k)),0) 
	dphi = lambda alpha: np.dot(np.transpose(p_k), (rosenbrock((x_k + alpha * p_k),1))) #derivative of phi 

	alpha_iminusone = 0
	alpha_i = 1 #first guess 
	alpha_max = 5
	i = 1 
	c1 = 1E-4
	c2 = 0.9

	#					Start Zoom function
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def zoom(alpha_low, alpha_high): 
		'''Interpolation (aka averaing) between two choices of alpha'''

		zoomcount = 1

		while True:

			#
			d1 = dphi(alpha_low) + dphi(alpha_high) - 3 * ((phi(alpha_low) - phi(alpha_high))/(alpha_low - alpha_high))
			d2 = np.sign(alpha_high - alpha_low) * np.sqrt(d1**2 - dphi(alpha_low) * dphi(alpha_high))

			alpha_j = alpha_high  - (alpha_high - alpha_low)* ((dphi(alpha_high) + d2 - d1)/(dphi(alpha_high) - dphi(alpha_low) + 2 * d2))
			#alpha_j = (alpha_high + alpha_low)/2.0 #simple midpoint 

			if (phi(alpha_j) > phi(0) + c1*alpha_j*dphi(0)) or (phi(alpha_j) > phi(alpha_low)): 
				alpha_high = alpha_j
				zoomcount += 1 
			else:
				if abs(dphi(alpha_j)) < -c2*dphi(0):
					alpha_star = alpha_j
					zoomcount += 1 
					break
				if dphi(alpha_j)*(alpha_high - alpha_low) >= 0:
					alpha_high = alpha_low
				alpha_low = alpha_j
				zoomcount += 1 

		return alpha_star,zoomcount


	#					End zoom function
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	#					Start finding alpha
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	while True: 

		if (phi(alpha_i) > phi(0) + c1*alpha_i*dphi(0)) or (phi(alpha_i) >= phi(alpha_iminusone) and i >1):
			alpha_star,zoomcount = zoom(alpha_iminusone, alpha_i)
			inner_counter += zoomcount
			break

		if abs(dphi(alpha_i)) <= -c2 *dphi(0): 
			alpha_star = alpha_i
			inner_counter += 1 
			break

		if dphi(alpha_i) >= 0:
			alpha_star,zoomcount = zoom(alpha_i,alpha_iminusone)
			inner_counter += zoomcount
			break

		alpha_iminusone = alpha_i
		alpha_i = alpha_i + 0.5


	#					End finding alpha
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	return alpha_star,inner_counter

#						End Definition for Strong Wolfe Search Algorithm 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#						Start Definition for BFGS function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def BFGS(algo):

	'''BFGS Pseudo Hessian that can use one of two linesearch methods: backtracking or H2 syle Strong Wolfe'''

	x_0 = rosenbrock(0,-1) #first argument doesn't matter, -1 return the initial conditions 
	H_0 = np.identity(len(x_0)) #Guessing identity matrix as inital guess for H
		
	epsilon = 1E-8 #threshold for stopping loop

	x_k = x_0 #setting current values to inital conditions 
	H_k = H_0 

	grad = rosenbrock(x_k,1) #Using provided (translated) rosenbrock funtion for evaluating gradient
	grad_norm = np.linalg.norm(grad)

	I = np.identity(len(x_k)) #Identity matrix I, needed for update to H_k+1

	#					Book Keeping 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	alpha_list = [] 
	function_values = []
	grad_norm_list = []
	criteria_list = []
	inner_counter = 0
	outer_counter = 1.0
	function_values.append(rosenbrock(x_k,0))
	grad_norm_list.append(grad_norm)
	time_start = time.time()
	#					End Book Keeping 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	#					Start BFGS 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	while grad_norm > epsilon: 

		p_k = np.dot(-H_k,grad)

		if algo == 'backtracking':
			alpha_k,temp_inner = backtracking_search(p_k,x_k,grad) #return step length, and inner counter
		elif algo == 'LS_Strong_Wolfe':
			alpha_k,temp_inner = LS_Strong_Wolfe(p_k,x_k,grad)

		print("Alpha = {}".format(alpha_k))


		x_kplusone = x_k + alpha_k * p_k #take step in p_k direction 

		s_k = x_kplusone - x_k

		y_k = rosenbrock(x_kplusone, 1) - grad

		rho_k = 1.0/np.dot(np.transpose(y_k),s_k)

		H_kplusone = (I - rho_k * np.dot(s_k,np.transpose(y_k))) * H_k * (I - rho_k * np.dot(y_k,np.transpose(s_k))) + rho_k * np.dot(s_k,np.transpose(s_k)) #BFGS update

		H_k = H_kplusone 

		x_k = x_kplusone 

		grad = rosenbrock(x_k,1) #update gradient and its norm 

		grad_norm = np.linalg.norm(grad)


						#Convergence Criteria (might not work for backtracking)
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		criteria_list.append((np.linalg.norm((H_k - rosenbrock(x_k,2)) * p_k)) / np.linalg.norm(p_k))
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


						#book keeping 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		grad_norm_list.append(grad_norm) #book keeping
		function_values.append(rosenbrock(x_k,0))
		outer_counter += 1
		inner_counter += temp_inner
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		assert np.dot(np.transpose(s_k),y_k) > 0, "Curvature condition failed!"

						#End BFGS
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	time_stop = time.time()
	function_values = [i[0][0] for i in function_values] #extract floats from a list of arrays, needed for plotting

	print("Done in {} seconds".format((time_stop - time_start)))


	plt.semilogy(function_values, label = "Objective value")
	plt.semilogy(grad_norm_list, label = "Norm of the gradient")
	plt.grid()
	plt.legend(loc = 'upper right')
	plt.xlabel("Iteration number")
	plt.title('BFGS {} Linesearch Results \n Outer Iterations: {} Inner Iterations: {}'.format(algo,outer_counter, inner_counter))
	plt.show()

#						End Definition for BFGS function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#						Start Full Newton Step Function Definition 
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

		print("Alpha = {}".format(alpha_k))

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
	#supply either 'backtracking' or 'LS_Strong_Wolfe' as an argument to either function. 

	BFGS('backtracking')
	#Full_Newton()