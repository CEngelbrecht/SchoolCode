import numpy as np 
import matplotlib.pyplot as plt 
import sympy as sp 
import math 

x1,x2,f = sp.symbols('x1,x2,f')

func = 100* (x2 - x1**2)**2 + (1 - x1)**2#Symbolic representation of f(x1,x2)

f = sp.lambdify((x1,x2),func,"numpy") #A lambdified version: speeds it up, faster than evalf(). Evaluate this with f(x1,x2)

grad_f = [sp.diff(func,x1),sp.diff(func,x2)] #Symbolic representation of the gradient 

grad_f = [sp.lambdify((x1,x2),grad_f[0],"numpy"),sp.lambdify((x1,x2),grad_f[1],"numpy")] #A list of lambdified functions

hessian = np.array([[sp.diff(sp.diff(func,x1),x1),sp.diff(sp.diff(func,x1),x2)], \
					[sp.diff(sp.diff(func,x2),x1),sp.diff(sp.diff(func,x2),x2)]])

hessian = np.array([[sp.lambdify((x1,x2),hessian[0][0],"numpy"),sp.lambdify((x1,x2),hessian[0][1],"numpy")],\
				    [sp.lambdify((x1,x2),hessian[1][0],"numpy"),sp.lambdify((x1,x2),hessian[1][1],"numpy")]])

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

def eval_hessian(x_k):

	p1 = x_k[0]
	p2 = x_k[1]

	h1 = hessian[0][0](p1,p2)
	h2 = hessian[0][1](p1,p2)
	h3 = hessian[1][0](p1,p2)
	h4 = hessian[1][1](p1,p2)

	return np.array([[h1,h2],[h3,h4]])


x_k_list = []
f_list = [] 

x_k = (1.2,1.2)
function_value = eval_function(x_k)

f_list.append(function_value)

l = 0 

while function_value > 1E-8:

	gradient = eval_grad(x_k)

	grad_norm = np.linalg.norm(eval_grad(x_k)) #returns the norm of the gradient at x_k 

	p_k = (-1.0) * (gradient / grad_norm) 

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Begin picking step length~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	
	def zoom(alpha_low, alpha_high):

		'''Zoom function that exists per iteration'''


		while True: 

		# 	#Try quadratic interpolation first

		 	#alpha_quad = -alpha_high**2 * dphi_0/ 2*(phi(alpha_high) - phi(alpha_low) - alpha_high*dphi_0) 
		# 	print("alpha_low = {}, alpha_high = {}, alpha_quad = {}".format(alpha_low,alpha_high,alpha_quad))


		# 	if phi(alpha_quad) <= phi_0 + c1 * alpha_quad * dphi_0: #Armijo Condition is unsatisfied from quad interp, d0 cubic

		# 		#cubic interpolation. Right now the subscripts are still alpha_1 and alpha_2, will change once I know quad works
		# 		print("Quad didn't work, doing cubic")

		# 		denom = alpha_0**2 * alpha_1**2 * (alpha_1 - alpha_0)
		# 		ar1 = np.array([[alpha_0**2, - alpha_1**2],[-alpha_0**3, alpha_1**3]]) #2x2 array
		# 		ar2 = np.array([[phi(alpha_1) - phi(0) - alpha_1 * dphi(0)],[phi(alpha_0) - phi(0) - alpha_0 * dphi(0)]]) #2x1 vector 

		# 		a,b = (1 / denom )*  np.dot(ar1,ar2) #solve linear system, get a and b coefficients 

		# 		alpha_cube = (-b + np.sqrt(b**2 - 3 * a * dphi(0)))/(3*a) #get cubic result

			alpha_interp =  ((alpha_low + alpha_high)/2.0)

			alpha_j = alpha_interp #Finished interpolating, set this to quad or cubic

			if (phi(alpha_j) > phi_0 + c1*alpha_j*dphi_0) or (phi(alpha_j) >= phi(alpha_low)): #line 4 of zoom

				alpha_high = alpha_j

			else: 

				if abs(dphi(alpha_j)) <= -c2*dphi_0: 

					alpha_star = alpha_j

		 			return alpha_star

		 		if dphi(alpha_j)*(alpha_high - alpha_low) >= 0: 

					alpha_high = alpha_low

				alpha_low = alpha_j
	
	phi_list = [] #book keeping 
	alpha_list = []

	alpha_0 = 0
	alpha_list.append(alpha_0)
	alpha_1 = 1
	alpha_max = 5 

	c1 = 1E-4
	c2 = 0.9 
	
	i = 1

	phi = lambda alpha: f(x_k[0] + alpha*p_k[0],x_k[1] + alpha*p_k[1]) #phi as a function of alpha, at the current x_k and p_k. Eval with phi(alpha)

	dphi = lambda alpha: np.dot(np.transpose(p_k), eval_grad(x_k + alpha * p_k)) #The derivative of phi. Eval with dphi(alpha)

	dphi_0 = dphi(0)

	phi_0 = phi(0)

	alpha = alpha_1

	while True:
			
		#print("i = {}".format(i))

		alpha_list.append(alpha) #book keeping

		phi_i = phi(alpha) #phi_i = phi(alpha_i)

		phi_list.append(phi_i)

		if (phi_i > phi_0 + c1*alpha*dphi_0) or (phi_i >= phi_list[i-1] and i > 0):

				print("phi_i = {}, phi_0 = {}, alpha = {}, dphi_0 = {}".format(phi_i,phi_0,alpha,dphi_0))

				print("alpha_list[i-1] = {}, alpha = {}".format(alpha_list[i-1],alpha))

				#zoom(alpha_low,alpha_high)				
				alpha_star = zoom(alpha_list[i-1],alpha)
				
				print "Breaking with alpha_star = {}".format(alpha_star)
				
				break #should break out of the while True loop

		if abs(dphi(alpha)) <= -c2*dphi_0: #line 7 of algorithm

				alpha_star = alpha

				break #should break out of the while True loop

		if (dphi(alpha) >= 0 and i > 0): #line 9 of algorithm

				alpha_star = zoom(alpha,alpha_list[i-1])

				break #should break out of the while True loop

		print('Increasing alpha')
		alpha = abs(alpha_max - alpha)*0.2


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

	l+= 1 
