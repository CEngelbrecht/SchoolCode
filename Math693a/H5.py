import numpy as np 
from rosenbrock_2Nd_translated import rosenbrock

alpha_list = [] 
function_value = []
inner_counter = 0 
outer_counter = 0 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def linesearch(p_k,x_k, grad): 

	alpha_bar = 1
	rho = 0.5
	c = 1E-4

	alpha = alpha_bar

	test_val1 = rosenbrock((x_k + (alpha * p_k)),0) #function evaluated at x_k + alpha * p_k 

	test_val2 = rosenbrock(x_k,0) + (c * alpha * np.transpose(p_k).dot(grad))

	while test_val1 > test_val2: 
		
		alpha = rho * alpha

		test_val1 = rosenbrock((x_k + (alpha * p_k)),0)

		test_val2 = rosenbrock(x_k,0) + (c * alpha * np.transpose(p_k).dot(grad))

		print test_val2,test_val1, alpha

	return alpha

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

x_0 = rosenbrock(0,-1) #first argument doesn't matter, -1 return the initial conditions 

H_0 = np.identity(len(x_0))
	
epsilon = 1E-8 #threshold 

x_k = x_0 #initial conditions 
H_k = H_0 

grad = rosenbrock(x_k,1)

I = np.identity(len(x_k)) #Initial guess to the inverse Hessian

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

while np.linalg.norm(grad) > epsilon: 

	p_k = np.dot(-H_k,grad)

	alpha_k = linesearch(p_k,x_k,grad) #return step length

	x_kplusone = x_k + alpha_k * p_k #take step in p_k direction 

	s_k = x_kplusone - x_k

	y_k = rosenbrock(x_kplusone, 1) - grad

	rho_k = 1.0/np.dot(np.transpose(y_k),s_k)

	H_kplusone = (I - rho_k * np.dot(s_k,np.transpose(y_k))) * H_k * (I - rho_k * np.dot(y_k,np.transpose(s_k))) + rho_k * np.dot(s_k,np.transpose(s_k)) #BFGS update

	H_k = H_kplusone

	x_k = x_kplusone

	grad = rosenbrock(x_k,1)



print "Done!"

print x_k