import matplotlib.pyplot as plt 
import numpy as np 
import sympy as sp 

x1,x2 = sp.symbols('x1 x2')

func = 10*(x2 - x1**2)**2 + (1 - x1)**2 

x_0 = (0,-1)

f = sp.lambdify((x1,x2),func,"numpy") #A lambdified version: speeds it up, faster than evalf(). Evaluate this with f(x1,x2)

hessian = np.array([[sp.diff(sp.diff(func,x1),x1),sp.diff(sp.diff(func,x1),x2)], \
					[sp.diff(sp.diff(func,x2),x1),sp.diff(sp.diff(func,x2),x2)]])

hessian = np.array([[sp.lambdify((x1,x2),hessian[0][0],"numpy"),sp.lambdify((x1,x2),hessian[0][1],"numpy")],\
				    [sp.lambdify((x1,x2),hessian[1][0],"numpy"),sp.lambdify((x1,x2),hessian[1][1],"numpy")]])

def eval_hessian(x_k):

	p1 = x_k[0]
	p2 = x_k[1]

	h1 = hessian[0][0](p1,p2)
	h2 = hessian[0][1](p1,p2)
	h3 = hessian[1][0](p1,p2)
	h4 = hessian[1][1](p1,p2)

	return np.array([[h1,h2],[h3,h4]])


