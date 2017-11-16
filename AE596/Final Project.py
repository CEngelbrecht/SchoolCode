import numpy as np 

V_xstar = 1594.6
V_ystar = 0.0 
m_0 = 17513 
g_0 = 9.81 
g_moon = 1.619
m_dot = -19.118 
Thrust = 60051.0
I_sp = - Thrust/(m_dot * g_0)
epsilon = 10E-2

def get_tgo(V_x,V_y,t_go):

	while 1: 

		dv = np.sqrt((V_xstar - V_x)**2 + (V_ystar - V_y + g_moon * t_go)**2)

		t_go_kplusone = (-m_0/abs(m_dot)) * (np.exp(-dv/(I_sp * g_0)) - 1)

		error = abs(t_go_kplusone - t_go)

		if error > epsilon:

			print("trying again with error {}".format(error))

			t_go += 1 

		else: 

			print("gotcha bitch")
			break 

	return t_go_kplusone

result = get_tgo(0, 0, 1)
print(result)