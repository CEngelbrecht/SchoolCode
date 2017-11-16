import numpy as np 
from rosenbrock_2Nd_translated import rosenbrock

x_0 = rosenbrock(0,-1) #first argument doesn't matter, -1 return the initial conditions 

H_0 = np.identity(len(x_0))
	
epsilon = 1E-8

x_k = x_0
H_k = H_0

grad = rosenbrock(x_k,1)

I = np.identity(len(x_k))

#while np.linalg.norm(grad) > 1E-8: 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def linesearch(p_k,x_k): 

	alpha_bar = 1
	rho = 0.5
	c = 1E-4

	alpha = alpha_bar

	test_val1 = rosenbrock((x_k + (alpha * p_k)),0)
	print test_val1

	while rosenbrock((x_k + alpha * p_k),0) > rosenbrock(x_k,0) + c * alpha * np.dot(np.transpose(p_k),rosenbrock(x_k,1)): 

		alpha = rho * alpha

	return alpha

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for i in range(1):

	p_k = np.dot(H_k,grad)

	alpha_k = linesearch(p_k,x_k)

	x_kplusone = x_k + alpha_k * p_k

	s_k = x_kplusone - x_k

	y_k = rosenbrock(x_kplusone, 1) - grad

	rho_k = 1.0/np.dot(np.transpose(y_k),s_k)

	H_kplusone = (I - rho_k * np.dot(s_k,np.transpose(y_k))) * H_k * (I - rho_k * np.dot(y_k,np.transpose(s_k))) + rho_k * np.dot(s_k,np.transpose(s_k))

	H_k = H_kplusone

	x_k = x_kplusone

	grad = rosenbrock(x_k,1)

