	d import numpy as np 

def rosenbrock(x, order): 

	if( order == - 1 ):
		'''Return initial conditions'''

		initial = [1.2000,1.2000,1.1000,1.1000,1.0500,1.0500,1.0250,1.0250,-1.2000,1.0000,-0.1000,1.0000,0.4500,1.0000,0.7250,1.0000,-2.4000,2.0000]

		R = np.zeros((len(initial),1))

		for k in range(len(R)):
			R[k] = initial[k]

		return R

	rb2d       = lambda x: (100* (x[1] - x[0]**2)**2 + (1 - x[0])**2)
	rb2d_x     = lambda x: (-400 * (x[1] - x[0]**2) * (x[0]) - 2 + 2 * x[0])
	rb2d_xx    = lambda x: (1200 * x[0]**2 - 400*x[1] + 2 )
	rb2d_xy    = lambda x: (-400 * x[0])
	rb2d_y     = lambda x: (200 * x[1] - 200 * x[0]**2)
	rb2d_yy    = lambda x: np.array(200)
	rb2d_grad  = lambda x: np.array([[rb2d_x(x)][0],[rb2d_y(x)][0]]) 
	rb2d_hess  = lambda x: np.array([[rb2d_xx(x),rb2d_xy(x)],[rb2d_xy(x),rb2d_yy(x)]])


	nx = len(x)

	if order == 0: 

		R = np.zeros((1,1))

		for k in range(0,nx,2):

			R = R + rb2d(x[k:(k+2)])

		return R

	elif order == 1: 

		R = np.zeros((len(x),1))

		for k in range(0,nx, 2):

			R[k:(k+2)] = rb2d_grad(x[k:(k+2)])

		return R 

	elif order == 2: 

		R = np.zeros((len(x),len(x)))

		for k in range(0, nx , 2):

			R[k:(k+2),k:(k+2)] = rb2d_hess(x[k:(k+2)])

		return R

	else: 
		print "rosenbrock couldn't use order {}".format(order)


