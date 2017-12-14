import matplotlib.pyplot as plt
from numpy.linalg import norm 
from numpy import linspace,meshgrid

def PLOTTER(list,globalstrat,tag,**kwargs):

	if 'delta_list' in kwargs:
		delta_list = kwargs['delta_list']

	if tag == 'function':
		plt.figure()
		list = [float(list[i]) for i in range(len(list))] #extract floats from list of 1x1 arrays
		plt.semilogy(list)
		plt.xlabel('iteration number')
		plt.title("Objective value")
		plt.grid()		

	if tag == 'direction' and globalstrat == 1:
		'''This only works for 2D case'''
		plt.figure()
		x_list = [i[0] for i in list]
		y_list = [i[1] for i in list]
		plt.plot(x_list,y_list,color = 'red',linewidth = 2,marker = '^')
		plt.grid()
		plt.xlabel("x")
		plt.ylabel('y')
		plt.title("x y position of 2D Rosenbrock")

		#plt.figure()
		xlist = linspace(-1.5,1.5,100)
		ylist = linspace(-1,3,100)
		X, Y = meshgrid(xlist, ylist)
		Z = 100 * (Y - X**2)**2 + (1 - X)**2
		cp = plt.contourf(X, Y, Z)

	elif tag == 'direction' and globalstrat == 3: 
		plt.figure()
		x_list = [i[0] for i in list]
		y_list = [i[1] for i in list]
		plt.plot(x_list,y_list,color = 'red',linewidth = 2,marker = '^')
		plt.grid()
		plt.xlabel("x")
		plt.ylabel('y')
		plt.title("x y position of 2D Rosenbrock")

		xlist = linspace(-1.5,1.5,100)
		ylist = linspace(-1,3,100)
		X, Y = meshgrid(xlist, ylist)
		Z = 100 * (Y - X**2)**2 + (1 - X)**2
		cp = plt.contourf(X, Y, Z)

	if tag == 'gradient':

		grad_norm = [norm(i) for i in list]

		plt.figure()
		plt.semilogy(grad_norm)
		plt.xlabel('iteration')
		plt.ylabel('gradient norm')
		plt.title("Norm of the gradient")

	plt.show()