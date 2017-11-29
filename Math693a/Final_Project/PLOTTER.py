import matplotlib.pyplot as plt 

def PLOTTER(list,tag):
	if tag == 'function':
		plt.figure()
		list = [float(list[i]) for i in range(len(list))] #extract floats from list of 1x1 arrays
		plt.semilogy(list)
		plt.xlabel('iteration number')
		plt.ylabel('Objective value')
		plt.grid()
		

	if tag == 'directon':
		'''This only works for 2D case'''

		pass
	plt.show()