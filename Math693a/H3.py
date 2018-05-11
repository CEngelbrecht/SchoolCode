import matplotlib.pyplot as plt 
import numpy as np 
import sympy as sp 

x1,x2 = sp.symbols('x1 x2')

func = 10*(x2 - x1**2)**2 + (1 - x1)**2 

x_0 = (0,-1)

f = sp.lambdify((x1,x2),func,"numpy") #A lambdified version: speeds it up, faster than evalf(). Evaluate this with f(x1,x2)

grad = np.array([[sp.diff(func,x1)],\
				 [sp.diff(func,x2)]]) #take derivatives 

grad = np.array([[sp.lambdify((x1,x2),grad[0],"numpy")],\
				 [sp.lambdify((x1,x2),grad[1],"numpy")]]) #array of functions to evaluate gradient; use with eval_grad 

hessian = np.array([[sp.diff(sp.diff(func,x1),x1),sp.diff(sp.diff(func,x1),x2)], \
					[sp.diff(sp.diff(func,x2),x1),sp.diff(sp.diff(func,x2),x2)]])

hessian = np.array([[sp.lambdify((x1,x2),hessian[0][0],"numpy"),sp.lambdify((x1,x2),hessian[0][1],"numpy")],\
				    [sp.lambdify((x1,x2),hessian[1][0],"numpy"),sp.lambdify((x1,x2),hessian[1][1],"numpy")]])

def eval_grad(x_k):

	p1 = x_k[0]
	p2 = x_k[1]

	g1 = grad[0][0](p1,p2)
	g2 = grad[1][0](p1,p2)

	return np.array((g1,g2)) #returns 2x1 array of evaluated gradient

def eval_hessian(x_k):

	p1 = x_k[0]
	p2 = x_k[1]

	h1 = hessian[0][0](p1,p2)
	h2 = hessian[0][1](p1,p2)
	h3 = hessian[1][0](p1,p2)
	h4 = hessian[1][1](p1,p2)

	return np.array([[h1,h2],[h3,h4]])

def full_step_pk(x_k): 

	return -1.0 * np.dot(np.linalg.inv(eval_hessian(x_k)) , eval_grad(x_k))

def steepest_descent_pk(x_k):

	grd = eval_grad(x_k)
	B_k = eval_hessian(x_k)

	return -1.0 * grd/np.linalg.norm(grd)

	#return (-1.0) * np.dot((np.dot((np.transpose(grd), grd)) / np.dot(np.dot(np.transpose(grd),B_k),grd)),grd) 

deltas = np.arange(0,2.2,0.2)  #values from 0 to 2 in steps of 0.2

#~~~~~~~~~~~~~~~Dogleg method~~~~~~~~~~~~~~~~~~~~~~~~~
dogleg_pk_list = [] 

for delta in deltas: 

	if np.linalg.norm(steepest_descent_pk(x_0)) >= delta: 

		dogleg_pk = delta * (steepest_descent_pk(x_0))/np.linalg.norm(steepest_descent_pk(x_0))

		dogleg_pk_list.append(dogleg_pk)

	elif np.linalg.norm(full_step_pk(x_0)) <= delta: 

		dogleg_pk = full_step_pk(x_0)

		dogleg_pk_list.append(dogleg_pk)
	else: 

		taus = np.arange(1,2,0.001)

		while 1: 
			print("trying to find tau")
			for tau in taus:

				if np.linalg.norm(steepest_descent_pk(x_0) + (tau - 1) * (full_step_pk(x_0) - steepest_descent_pk(x_0)))**2 == delta:

					break 
					print("found tau")
				else:
					 pass

		dogleg_pk = steepest_descent_pk + (tau - 1) * (full_step_pk(x_0) - steepest_descent_pk(x_0))

		dogleg_pk_list.append(dogleg_pk)


print len(dogleg_pk_list)
print len(deltas)
#Exact trust region algorithm: 


