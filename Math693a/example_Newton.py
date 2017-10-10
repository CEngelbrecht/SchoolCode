import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 
import sympy as sp 

x_0 = [1.0,1.0]

x1,x2,f = sp.symbols('x1,x2,f')

f = (x1 + x2**2)**2 #Symbolic representation of f(x1,x2)

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
	p2_2 = hessian[0][1].subs([(x1,x_k[0]),(x2,x_k[1])]).evalf()

	return np.array([[p1_1,p1_2],[p2_1,p2_2]])


x_k = x_0 

x_k_list = []
#x_k_list.append(x_k)

function_value = eval_function(x_0) #inital value of the function

f_list = []
#f_list.append(function_value)

alpha_list = []

while function_value > 1E-8:
	
	gradient = eval_grad(x_k)

	evaluated_hessian = eval_hessian(x_k)

	p_k = (-1) * np.linalg.inv(evaluated_hessian.astype(np.float64)).dot(gradient) #Take the inverse of the evaluated hessian typecasted to float
																				   # and take the dot product of that with the gradient

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




#Plot original function 

fig = plt.figure()
ax = fig.gca(projection='3d')

X = np.arange(-1.5, 1.2, 0.01)
Y = np.arange(-1.5, 1.2, 0.01)
X, Y = np.meshgrid(X, Y)

Z = (X + Y**2)**2 

surf = ax.plot_surface(X, Y, Z,cmap=cm.coolwarm,linewidth=-1	, antialiased=True)

#plot the x_k values 

for i in range(len(x_k_list)): 
	ax.scatter(x_k_list[i][0],x_k_list[i][1],f_list[i],c= 'black',marker = '^',s = 100)# place down x_k markers

xs = [x_k_list[i][0] for i in range(len(x_k_list))]
ys = [x_k_list[i][1] for i in range(len(x_k_list))]

ax.plot(xs,ys,f_list,c = 'red')

plt.show()