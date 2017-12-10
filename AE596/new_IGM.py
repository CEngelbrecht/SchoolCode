def bashars_IGM(y,V_x,V_y,m,t_go):

	from numpy import sin,cos,log,array,arctan
	from numpy.linalg import solve 

	m_dot = -19.118
	y_f = 185200
	g_0 = 9.81
	g_m = 1.619
	T = 60051.0
	I_sp = - T/(m_dot * g_0)
	V_ystar = 0
	V_xstar = 1594.6

 
	chi_tilde = arctan(((V_ystar - V_y) - (g_m * t_go) )/ ((V_xstar - V_x)))

	S = (m + m_dot * t_go)/m

	b1 = y_f - y - (V_y * t_go) + g_m/2.0 *t_go - ((T*sin(chi_tilde) * m)/abs(m_dot)**2) * (S*log(abs(S)) - S + 1.0)

	a_11 = ((T * cos(chi_tilde)*m**2)/abs(m_dot)**3) * (S * log(abs(S)) - S**2/2.0 +0.5)

	a_12 = (T * cos(chi_tilde) * m)/abs(m_dot)**2 * (S * log(abs(S)) - S + 1)

	b2 = -I_sp * g_0 *log(S) * sin(chi_tilde) + (T*sin(chi_tilde)/abs(m_dot)) * log(abs(S))

	a_21 = -(T*cos(chi_tilde)/m_dot**2) * m *log(abs(S)) + abs(m_dot)*t_go

	a_22 = -(T*cos(chi_tilde)/abs(m_dot)) * log(abs(S))

	a_array = array([[a_11,a_12],[a_21,a_22]])

	b_array = array([[b1],[b2]])

	k1,k2 = solve(a_array,b_array) 

	return chi_tilde + k2

