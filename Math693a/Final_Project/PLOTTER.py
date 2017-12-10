import matplotlib.pyplot as plt
from numpy.linalg import norm 

def PLOTTER(list,tag):

	if tag == 'function':
		plt.figure()
		list = [float(list[i]) for i in range(len(list))] #extract floats from list of 1x1 arrays
		plt.semilogy(list)
		plt.xlabel('iteration number')
		plt.ylabel('Objective value')
		plt.grid()		

	if tag == 'direction':
		'''This only works for 2D case'''
		plt.figure()
		x_list = [i[0] for i in list]
		y_list = [i[1] for i in list]
		plt.plot(x_list,y_list)
		plt.grid()
		plt.xlabel("x")
		plt.ylabel('y')
		plt.title("x y position of 2D Rosenbrock")

	if tag == 'gradient':

		grad_norm = [norm(i) for i in list]

		plt.figure()
		plt.semilogy(grad_norm)
		plt.xlabel('iteration')
		plt.ylabel('gradient norm')
		plt.title("Norm of the gradient")

	plt.show()