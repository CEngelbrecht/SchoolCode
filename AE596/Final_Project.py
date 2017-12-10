'''
Christian Engelbrecht
AE596 Final Project 
December 7th 2017
'''

from numpy.linalg import solve
from numpy import sin,cos,log,array,sqrt,exp,arctan
import time 
import matplotlib.pyplot as plt

def get_tgo(V_x,V_y,m):

	t_guess = 50
	t_go = t_guess
	epsilon = 1E-4
	V_xstar = 1700 #m/s
	V_ystar = 0
	g_m = 1.619
	g_0 = 9.81 #m/2 
	m_dot = -19.118
	T = 60051.0
	I_sp = - T/(m_dot * g_0)

	while True:

		dv = sqrt((V_xstar - V_x)**2 + (V_ystar - V_y + g_m * t_go)**2)

		t_go_kplusone = (-1.0*m/abs(m_dot)) * (exp(-dv/(I_sp * g_0)) - 1)

		error = abs(t_go - t_go_kplusone)

		if error > epsilon:

			#print("trying again with error {}".format(error))

			t_go = t_go_kplusone

			if error > 1E3: 

				print("Not converging")
				break

		else: 

			break 

	return t_go_kplusone

def IGM_GUIDANCE(y,V_x,V_y,m,t_go): 

	chi_tilde = arctan(((V_ystar - V_y) - (g_m * t_go) )/ ((V_xstar-V_x)))

	a_20 = V_y +((T * sin(chi_tilde))/abs(m_dot)) * log(m / (m - (abs(m_dot)*t_go))) - (g_m * t_go)

	a_21 = T * cos(chi_tilde) * ( (m * log((m)/(m - abs(m_dot)*(t_go))))/m_dot**2 - t_go/abs(m_dot))

	a_22 = (T * cos(chi_tilde) / (abs(m_dot))) * log(m/(m + (abs(m_dot)*t_go)))

	a_10 = y + (V_y * t_go) - (g_m * (t_go**2)/2.0)  + ( (T * sin(chi_tilde))/(abs(m_dot))) *(t_go + ((t_go -(m/abs(m_dot)))) * log(m/(m - (abs(m_dot) * t_go))))

	a_11 = (T * cos(chi_tilde)/(abs(m_dot)**2)) * ( (-(abs(m_dot) * t_go**2)/ 2.0) + (m * ( (t_go - (m/abs(m_dot)))* log(m/(m - (abs(m_dot)*t_go))) + t_go)))

	a_12 = (T * cos(chi_tilde)/ (abs(m_dot))) * ( (t_go - (m/abs(m_dot))) * log(m / (m - (abs(m_dot) * t_go))) + t_go);

	b2 = I_sp*g_0 * log((m)/(m - abs(m_dot) * t_go)) * sin(chi_tilde) + V_y - g_m * t_go - a_20
	b1 = y_fstar + V_ystar * t_go - a_10

	a_array = array([[a_11,a_12],[a_21,a_22]])
	b_array = array([[b1],[b2]])

	k1,k2 = solve(a_array,b_array)

	beta = chi_tilde + k2 

	return beta


def integrate_EOMS(x,y,V_x,V_y,m,beta):

	a = T/m

	x += V_x * dt 

	y += V_y * dt 

	V_x += a*cos(beta) * dt 

	V_y += (a*sin(beta) - g_m)  * dt 

	m = m - abs(m_dot) * dt

	return x,y,V_x,V_y,m
if __name__ == '__main__': 
	
	#~~~~~~~~~~~~~~~~~~
	x_0 = 0
	y_0 = 0 
	V_x_0 = 0 
	V_y_0 = 0 
	m_0 = 17513
	m_dot = -19.118 
	T = 60051.0
	g_m = 1.619
	g_0 = 9.81 #m/2
	I_sp = - T/(m_dot * g_0)
	dt = 1
	t = 0
	V_ystar = 0 

	V_xstar = 1700 #m/s
	y_fstar = 100000 #m 
	#~~~~~~~~~~~~~~~~~~

	x = x_0
	y = y_0
	V_x = V_x_0
	V_y = V_y_0
	m = m_0

	#~~~~~~~~~~~~~~~~~~Bookkeeping~~~~~~~~~~~~~
	x_list = []
	y_list = []
	V_x_list = [] 
	V_y_list = [] 
	m_list = []
	beta_list = []
	t_go_list = []
	t_list = []
	#~~~~~~~~~~~~~~~~~~Bookkeeping~~~~~~~~~~~~~

	t_go = get_tgo(V_x,V_y,m)
	t_go_list.append(t_go)

	while dt <= t_go:

		t_go = get_tgo(V_x,V_y,m)

		beta = float(IGM_GUIDANCE(y,V_x,V_y,m,t_go))

		x, y, V_x, V_y,m  = integrate_EOMS(x,y,V_x,V_y,m,beta)

		print x, y, V_x, V_y,m,beta

		t = t + dt 

		#~~~~~~~~~~~~~~~~~~Bookkeeping~~~~~~~~~~~~~
		x_list.append(x)
		y_list.append(y)
		V_x_list.append(V_x)
		V_y_list.append(V_y)
		m_list.append(m)
		beta_list.append(beta)
		t_go_list.append(t_go)
		t_list.append(t)
		#~~~~~~~~~~~~~~~~~~Bookkeeping~~~~~~~~~~~~~


	#final time 

	x, y, V_x, V_y,m  = integrate_EOMS(x,y,V_x,V_y,m,beta)

	#Plotting routines
	directory = "/home/rp/Documents/Christian_School/Plots/" 

	plt.plot(x_list,y_list)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.grid()
	plt.title("Mission 2 Trajectory")
	plt.savefig(directory+"M2_XvY")

	plt.figure()
	plt.plot(t_list,beta_list)
	plt.title("Mission 2" + r"$\beta$(t)")
	plt.ylabel('radian')
	plt.grid()
	plt.savefig(directory+"M2_BetavT")

	plt.figure()
	plt.title("Mission 2 y(t)")
	plt.plot(t_list,y_list)
	plt.grid()
	plt.savefig(directory+"M2_Yvt")
	
	plt.figure()
	plt.plot(t_list,[sqrt(i**2 + j**2) for i,j in zip(V_x_list,V_y_list)])
	plt.title("Mission 2 Velocty v Time")
	plt.ylabel("m/s")
	plt.grid()
	plt.savefig(directory+"M2_Vvt")

	plt.figure()
	plt.title("Mission 2" + r"$\gamma$(t)")
	plt.grid()
	plt.plot(t_list,[arctan(i/j) for i,j in zip(V_x_list,V_y_list)])
	plt.ylabel('radian')
	plt.savefig(directory+"M2_gammavt")

	plt.show()



