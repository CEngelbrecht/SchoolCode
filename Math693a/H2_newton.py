import numpy as np 
import matplotlib.pyplot as plt 
import sympy as sp 
import math 
import time 

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
alpha_list_outer = []

loc1 = (1.2,1.2)
loc2 = (-1.2,1.0)
loc = loc2
#loc = loc1

#Choosing location starting point 

if loc == loc1:
	x_0 = loc1
elif loc == loc2:
	x_0 = loc2
x_k = x_0

function_value = eval_function(x_k)

gradient = eval_grad(x_k)

f_list.append(function_value)


while np.linalg.norm(gradient) > 1E-12:

	gradient = eval_grad(x_k)

	evaluated_hessian = eval_hessian(x_k)

	strat = "SD"
	#strat = "Newton"

	if strat == "Newton":
		p_k = (-1) * (np.linalg.inv(evaluated_hessian.astype(np.float64))).dot(gradient) #Take the inverse of the evaluated hessian typecasted to float
	elif strat == "SD":
		grad_norm = np.linalg.norm(gradient) #returns the norm of the gradient at x_k 
		p_k = (-1.0) * (gradient / grad_norm) 

	inner_counter = 0 

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Begin picking step length~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

	#~~~~~~~~~~~~~~~Begin define zoom~~~~~~~~~~~~~~~~	
	def zoom(alpha_low, alpha_high):
		alpha_list_zoom = []
		j = 0 
		e_1 = 1E-3
		e_2 = e_1

		counter = 0 
		'''Zoom function that exists per iteration'''
		
		while True:
			counter += 1 
			

			#interpolate 

			d1 = dphi(alpha_low) + dphi(alpha_high) - 3 * ((phi(alpha_low) - phi(alpha_high))/(alpha_low - alpha_high))
			d2 = np.sign(alpha_high - alpha_low) * np.sqrt(d1**2 - dphi(alpha_low) * dphi(alpha_high))

			alpha_j = alpha_high  - (alpha_high - alpha_low)* ((dphi(alpha_high) + d2 - d1)/(dphi(alpha_high) - dphi(alpha_low) + 2 * d2))
			
			alpha_list_zoom.append(alpha_j)	

			if (phi(alpha_j) > phi(0) + c1*alpha_j*dphi_0) or (phi(alpha_j) >= phi(alpha_low)):

				alpha_high = alpha_j
			else: 
				if abs(dphi(alpha_j)) <= -c2*dphi(0):
					alpha_star = alpha_j
					break
				if dphi(alpha_j)*(alpha_high - alpha_low):
					alpha_high = alpha_low
				alpha_low = alpha_j 

			#safeguards
			if j > 0: 
	
				if abs(alpha_list_zoom[j] - alpha_list_zoom[j-1]) < e_1 or abs(alpha_list_zoom[j]) < e_2:
					print "Triggered a safeguard"
					alpha_star = alpha_list_zoom[j-1]/2.0  
					break
			j += 1

		return alpha_star,alpha_list_zoom
	#~~~~~~~~~~~~~~~End define zoom~~~~~~~~~~~~~~~~~~
	 
	alpha_list_inner = []
	alpha_0 = 0
	alpha_list_inner.append(alpha_0)
	alpha_1 = 1
	alpha_max = 5 

	c1 = 1E-4
	c2 = 0.9 
	
	i = 0

	phi = lambda alpha: f(x_k[0] + alpha*p_k[0],x_k[1] + alpha*p_k[1]) #phi as a function of alpha, at the current x_k and p_k. Eval with phi(alpha)

	dphi = lambda alpha: np.dot(p_k, eval_grad(x_k + alpha * p_k)) #The derivative of phi. Eval with dphi(alpha)

	dphi_0 = dphi(0)

	phi_0 = phi(0)

	alpha = alpha_1

	while True:

		inner_counter += 1 
			
		alpha_list_inner.append(alpha) #book keeping

		if i > 0:

			if (phi(alpha_list_inner[i]) > (phi(0) + c1*alpha_list_inner[i])) or (phi(alpha_list_inner[i]) > phi(alpha_list_inner[i-1])):
				alpha_star,alpha_list_zoom = zoom(alpha_list_inner[i-1],alpha_list_inner[i])
								
				break #should break out of the while True loop

		if abs(dphi(alpha_list_inner[i])) <= -c2*dphi(0): #line 7 of algorithm
			alpha_star = alpha_list_inner[i]

			break #should break out of the while True loop

		if (dphi(alpha_list_inner[i]) >= 0): #line 9 of algorithm
			alpha_star,alpha_list_zoom = zoom(alpha_list_inner[i],alpha_list_inner[i-1])

			break #should break out of the while True loop

		alpha += abs(alpha_max - alpha)*0.2  #increasing alpha
		i += 1 


	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End picking step length~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

	alpha_k = alpha_star #This is the alpha to use for this step in the process
	
	x_k = x_k + alpha_k*p_k

	alpha_list_outer.append(alpha_k)
	function_value = eval_function(x_k)
	gradient = eval_grad(x_k)
	x_k_list.append(x_k)
	f_list.append(function_value)
	

plt.semilogy(range(len(f_list)),f_list,'o-',c = 'black')
plt.grid()
plt.title(strat+" Strategy, H2 with Safeguards, x_0 = "+str(x_0))
plt.ylabel("Objective Value")
plt.show()