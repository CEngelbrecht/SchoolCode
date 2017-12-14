import numpy as np 


'''x0,x1,x2,x3 = sp.symbols("x0,x1,x2,x3")
	f = (x0 + 10*x1)**2 + 5*(x2 - x3)**2 + (x1 - 2*x2)**4 + 10*(x0 - x3)**4
	x_0 = (3,-1,0,1) #inital starting tuple
	x_0 = np.reshape(np.array(x_0),(len(x_0),1)) #numpy array of (x_0)T
	'''

def Powell_function(x,order):

	powell = lambda x: (x[0] + 10*x[1])**2 + 5*(x[2] - x[3])**2 + (x[1] - 2*x[2])**4 + 10*(x[0] - x[3])**4 
	
	powell_x0 = lambda x: 2*x[0] + 20*x[1] + 40*(x[0] - x[3])**3
	powell_x1 = lambda x: 20*x[0] + 200*x[1] + 4*(x[1] - 2*x[2])**3
	powell_x2 = lambda x: 10*x[2] - 10*x[3] - 8*(x[1] - 2*x[2])**3
	powell_x3 = lambda x: -10*x[2] + 10*x[3] - 40*(x[0] - x[3])**3

	powell_x0x0 = lambda x:  2*(60*(x[0] - x[3])**2 + 1)
	powell_x0x1 = lambda x:  20
	powell_x0x2 = lambda x:  0
	powell_x0x3 = lambda x:  -120*(x[0] - x[3])**2

	powell_x1x0 = lambda x:  20
	powell_x1x1 = lambda x:  4*(3*(x[1] - 2*x[2])**2 + 50)
	powell_x1x2 = lambda x:  -24*(x[1] - 2*x[2])**2
	powell_x1x3 = lambda x:  0  

	powell_x2x0 = lambda x:  0
	powell_x2x1 = lambda x:  -24*(x[1] - 2*x[2])**2
	powell_x2x2 = lambda x:  2*(24*(x[1] - 2*x[2])**2 + 5)
	powell_x2x3 = lambda x:  -10 

	powell_x3x0 = lambda x:  -120*(x[0] - x[3])**2
	powell_x3x1 = lambda x:  0 
	powell_x3x2 = lambda x:  -10
	powell_x3x3 = lambda x:  10*(12*(x[0] - x[3])**2 + 1)	

	powell_grad = lambda x: np.array([[powell_x0(x)[0]],[powell_x1(x)[0]],[powell_x2(x)[0]],[powell_x3(x)[0]]])

	powell_hess = lambda x: np.array([[powell_x0x0(x), powell_x0x1(x), powell_x0x2(x), powell_x0x3(x)], \
									  [powell_x1x0(x), powell_x1x1(x), powell_x1x2(x), powell_x1x3(x)], \
									  [powell_x2x0(x), powell_x2x1(x), powell_x2x2(x), powell_x2x3(x)], \
									  [powell_x3x0(x), powell_x3x1(x), powell_x3x2(x), powell_x3x3(x)]])

	nx = len(x)

	if order == 0: 

		R = np.zeros((1,1))

		for k in range(0,nx,4):

			R = R + powell(x[k:(k+4)])

		return R

	elif order == 1: 

		R = np.zeros((len(x),1))

		for k in range(0,nx,4):

			R[k:(k+4)] = powell_grad(x[k:(k+4)])

		return R

	elif order == 2: 

		R = np.zeros((len(x),len(x)))

		for k in range(0, nx , 4):

			R[k:(k+4),k:(k+4)] = powell_hess(x[k:(k+4)])

		return R