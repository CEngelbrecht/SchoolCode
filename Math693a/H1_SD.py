import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 
import sympy as sp 

x_0 = [-1.2,1.0]

x1,x2,f = sp.symbols('x1,x2,f')

f = 100* (x2 - x1**2)**2 + (1 - x1)**2#Symbolic representation of f(x1,x2)

grad_f = [sp.diff(f,x1),sp.diff(f,x2)] #Symbolic representation of the gradient 

hessian = np.array([[sp.diff(sp.diff(f,x1),x1),sp.diff(sp.diff(f,x1),x2)], \
					[sp.diff(sp.diff(f,x2),x1),sp.diff(sp.diff(f,x2),x2)]])


def eval_grad(x_k):
	'''Takes in two points, and evaluates the predefined symbolic represetation of the gradient at those points
		and returns the subsequent tuple'''

	p1 = grad_f[0].subs([(x1,x_k[0]),(x2,x_k[1])]).evalf()
	p2 = grad_f[1].subs([(x1,x_k[0]),(x2,x_k[1])]).evalf()

	return (p1,p2)

def eval_function(x_k):
	''' #Returns tuple after evaluating symbolic representation of f at x_k'''

	return f.subs([(x1,x_k[0]),(x2,x_k[1])]).evalf()

def eval_hessian(x_k):
	'''Returns 2x2 numpy array of the symbolic hessian function evaluated at (point1,point2)'''

	p1_1 = hessian[0][0].subs([(x1,x_k[0]),(x2,x_k[1])]).evalf()
	p1_2 = hessian[0][1].subs([(x1,x_k[0]),(x2,x_k[1])]).evalf()
	p2_1 = hessian[1][0].subs([(x1,x_k[0]),(x2,x_k[1])]).evalf()
	p2_2 = hessian[1][1].subs([(x1,x_k[0]),(x2,x_k[1])]).evalf()

	return np.array([[p1_1,p1_2],[p2_1,p2_2]])


x_k = x_0 

x_k_list = []
#x_k_list.append(x_k)

function_value = eval_function(x_0) #inital value of the function

f_list = []
#f_list.append(function_value)

#Steepest descent method: 
alpha_list = []
while function_value > 1E-8:
	
	gradient = eval_grad(x_k)

	grad_norm = np.sqrt(float(gradient[0]**2 + gradient[1]**2))

	p_k = (-1.0) * (gradient / grad_norm)

	#~~~~~~~~Backtracking~~~~~~~~~~~~~~~~~#

	alpha_bar = 1
	rho = 0.5
	c = 1E-4
	alpha = alpha_bar

	function_value = eval_function(x_k + p_k * alpha)

	test_value = eval_function(x_k) + c * alpha * p_k.dot(gradient)

	while function_value > test_value: 

		alpha = rho * alpha

		function_value = eval_function(x_k + p_k * alpha)

		test_value = eval_function(x_k) + c * alpha * p_k.dot(gradient)

	#~~~~~~~~End Backtracking~~~~~~~~~~~~~~~~~#


	alpha_k = alpha #This is the alpha to use for this step in the process

	function_value = eval_function(x_k) #regaining current value of the function for logging purposes

	print("alpha_k = {}, x_k = {}, function_value = {}, p_k = {}".format(alpha_k,x_k,function_value,p_k)) #Correct order of things 
	x_k_list.append(x_k)
	f_list.append(function_value)
	alpha_list.append(alpha)
	
	x_k = x_k + alpha_k*p_k

	x_k_list.append(x_k)
	f_list.append(function_value)


plt.subplot(121)
plt.semilogy(range(len(f_list)),f_list,'o-',c = 'black')
plt.grid()
plt.xlabel('Iteration number')
plt.ylabel('Objective Value')
plt.title('Function Value vs Iteration')

plt.subplot(122)
plt.scatter(range(len(alpha_list)),alpha_list)
plt.title('Step Size')
plt.xlabel('Iteration Number')
plt.ylabel('Step Size '+r'$\alpha$')

plt.suptitle('(-1.2,1.0) Rosenbrock Steepest Descent \n with Backtracking Line Search')

plt.show()