import numpy as np 
import matplotlib.pyplot as plt 
import sympy as sp 
import math 

x1,x2,f = sp.symbols('x1,x2,f')

f = 100* (x2 - x1**2)**2 + (1 - x1)**2#Symbolic representation of f(x1,x2)

grad_f = [sp.diff(f,x1),sp.diff(f,x2)] #Symbolic representation of the gradient 

f = sp.lambdify((x1,x2),f,"numpy") #A lambdified version: speeds it up, faster than evalf(). Evaluate this with f(x1,x2)

grad_f = [sp.lambdify((x1,x2),grad_f[0],"numpy"),sp.lambdify((x1,x2),grad_f[1],"numpy")] #A list of lambdified functions

def eval_grad(x_k):
	'''Returns a [list] of values where g1,g2 are the gradient components evaluated at (xk)'''

	p1 = x_k[0]
	p2 = x_k[1]

	g1 = grad_f[0](p1,p2)
	g2 = grad_f[1](p1,p2)

	return [g1,g2]

def eval_function(x_k):
	'''Returns scalar value of the funcion at x_k'''

	p1 = x_k[0]
	p2 = x_k[1]

	return f(p1,p2)


x_k_list = []
f_list = [] 

x_k = (1.2,1.2)
function_value = eval_function(x_k)

f_list.append(f)

while function_value > 1E-8:

	gradient = eval_grad(x_k)

	grad_norm = np.linalg.norm(eval_grad(x_k)) #returns the norm of the gradient at x_k 

	p_k = (-1.0) * (gradient / grad_norm) 

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Begin picking step length~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	
	def zoom(alpha_low, alpha_high):
		'''Zoom function that exists per iteration'''

		alpha_quad = -alpha_low**2 * dphi_0/ 2*(phi(alpha_low) - phi_0 - alpha_low*dphi_0) #quadratic interpolation

		if phi(alpha_quad) <= phi_0 + c1 * alpha_quad * dphi_0: #Armijo Condition is unsatisfied from quad interp, dp cubic

			#cubic interpolation, how do I do this? 
			pass

		alpha_j = alpha_quad #Finished interpolating, set this to quad or cubic

		if (phi(alpha_j) > phi_0 + c1*alpha_j*dphi_0) or (phi(alpha_j) >= phi(alpha_low)): #line 4 of zoom

			alpha_high = alpha_j

		else: 

			if dphi(alpha_j) <= -c2*dphi_0: 

				alpha_star = alpha_j

				return alpha_star

			if dphi(alpha_j)*(alpha_high * alpha_low) >= 0: 

				alpha_high = alpha_low

			alpha_high = alpha_low
	
	alpha_0 = 1
	c1 = 1E-4
	c2 = 0.9 
	
	i = 0 

	phi_list = []
	alpha_list = []

	phi = lambda alpha: f(x_k[0] + (alpha*p_k)[0],x_k[1] + (alpha*p_k)[1]) #phi as a function of alpha, at the current x_k and p_k. Eval with phi(alpha)

	dphi = lambda alpha: np.dot(p_k, eval_grad(x_k + alpha * p_k)) #The derivative of phi. Eval with dphi(alpha)

	dphi_0 = dphi(0)

	phi_0 = phi(0)

	alpha = alpha_0
	alpha_list.append(alpha)

	while True:
		print i 
		phi_i = phi(alpha)
		phi_list.append(phi_i)

		if i > 0: # just checking the first iteration 

			if (phi_i > dphi_0 + c1*alpha*dphi_0) or (phi_i >= phi_list[i-1]):
				print "one"

				alpha_star = zoom(alpha_list[i-1],alpha)
				print "Breaking with alpha_star = {}".format(alpha_star)
				break #should break out of the while True loop

			if abs(dphi(alpha)) <= -c2*dphi_0: #line 7 of algorithm
				print "two"

				alpha_star = alpha

				break #should break out of the while True loop

			if dphi(alpha) >= 0: 
				print "three"

				alpha_star = zoom(alpha,alpha_list[i-1])

				break #should break out of the while True loop
		else: 
			#Can't do lines 4 or 9, but only 7? 

			if abs(dphi(alpha)) <= -c2*dphi_0: #line 7 of algorithm

				alpha_star = alpha

				break 
		i += 1 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End picking step length~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

	alpha_k = alpha_star #This is the alpha to use for this step in the process

	function_value = eval_function(x_k) #regaining current value of the function for logging purposes

	print("alpha_k = {}, x_k = {}, function_value = {}, p_k = {}".format(alpha_k,x_k,function_value,p_k)) #Correct order of things 
	x_k_list.append(x_k)
	f_list.append(function_value)
	alpha_list.append(alpha)

	print alpha_k
	
	x_k = x_k + alpha_k*p_k

	x_k_list.append(x_k)
	f_list.append(function_value)
