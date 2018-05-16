from numpy import array,pi,deg2rad,sin,cos,dot,zeros,matmul
from numpy.linalg import inv,solve
import matplotlib.pyplot as plt 

def VOR2D(gamma,x,z,x_j,z_j):
	'''
	Calcualte velocity at x,z due to vortex of strength gamma located at x_j,z_j
	'''

	r_j = (x - x_j)**2 + (z - z_j)**2 
	u = gamma*(z - z_j) / (2*pi*r_j**2)
	w = -gamma*(x-x_j)/(2*pi*r_j**2)

	return u,w
	
def plot(x_colloc_locs,z_colloc_locs,vortex_x_locs,vortex_z_locs,TE_loc,LE_loc,*wake_loc):

	plt.plot(x_colloc_locs,z_colloc_locs,'g+',label = 'collation points')
	plt.plot(vortex_x_locs,vortex_z_locs,'ro',label = 'lumped vortices')
	plt.plot(TE_loc[0],TE_loc[1],'b^')
	plt.plot(LE_loc[0],LE_loc[1],'b^')
	plate_shape = [LE_loc[0],TE_loc[0]],[LE_loc[1],TE_loc[1]]
	if wake_loc:
		plt.plot(wake_loc[0][0],wake_loc[0][1],'r*')
	plt.plot(plate_shape[0],plate_shape[1],linewidth = 0.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

	#Initial Setup 
	c = 1.0
	dt = 0.1 
	T = 5 
	alpha = 5
	alpha = deg2rad(alpha)

	normal = array([sin(alpha),cos(alpha)])
	Q_inf = (0.05 * c)/dt 

	x_colloc_locs = [i + cos(alpha)*0.75*(c/4) for i in [0,0.25,0.5,0.75]]
	z_colloc_locs = [-i * sin(alpha) for i in x_colloc_locs]
	vortex_x_locs = [i + cos(alpha)*0.25*(c/4) for i in [0,0.25,0.5,0.75]]
	vortex_z_locs = [-i * sin(alpha) for i in vortex_x_locs]
	LE_loc = [0,0]
	TE_loc = [c*cos(alpha),-c*sin(alpha)]

	wake_locs = [TE_loc[0] + 0.25*Q_inf*dt*cos(alpha),TE_loc[1] - 0.25*Q_inf*dt*sin(alpha)]

	vortex_x_locs.append(wake_locs[0])
	vortex_z_locs.append(wake_locs[1])

	gamma_foil_list = []

	gamma = 1 

	A = zeros((5,5))
	B = zeros(5)

	# At t = dt, assume gamma = 1.0. 
	# Set up matrix for A * gamma = RHS d 
	# A matrix has influence coefficients a_11,a_12, etc 
	# Loop through all location points, call VOR2D at each point 

	for i in range(len(x_colloc_locs)):

		for j in range(len(x_colloc_locs)+1):

			u,w = VOR2D(gamma,x_colloc_locs[i],z_colloc_locs[i],vortex_x_locs[j],vortex_z_locs[j])

			A[i,j] = dot(array([u,w]), normal)


	A[4][:] = 1.0 
	gamma_foil0 = matmul(inv(A),B)
	gamma_sum = sum(gamma_foil0)
	gamma_foil_list.append(gamma_sum)

	# the above gamma is only for t = 0 
	# at t = dt, the plate is moved left at speed Q_inf, so new 
	# x_locs = old_x_loc - Q_inf*dt


		
	ITs = 0 #number of iterations 
	# wake_locs = zeros((ITs,2))
	# #wake_locs[0] = [TE_loc[0] + 0.25*Q_inf*dt*cos(alpha),TE_loc[1] - 0.25*Q_inf*dt*sin(alpha)] #first wake location, 0.2

	gamma_wake = 0
	gamma_wake_list = []
	gamma_wake_list.append(gamma_wake)
	wake_counter = 1 #introduce one wake vortex 

	for k in range(ITs):

		#loop over all collocation points 
		for i in range(4): 

			wake_velocies = zeros((wake_counter,2))
			RHS = zeros(5)
			vel_tot = zeros((n,2))

			for w in range(wake_counter):
			#loop as many times as there are wake vortices and determine induced velocity from all the wakes on the i'th collocation point 
				wake_velocites[w] = VOR2D(gamma_wake,x_colloc_locs[i],z_colloc_locs[i],vortex_x_locs[w+n],vortex_z_locs[w+n]) 

			vel_tot[i] = (sum(wake_velocity[:][0],sum(wake_velocity[:][1]))

			RHS[i] = dot((-Q_inf + vel_tot[i]),normal)

		for i in range(4):

			RHS[i] = dot([-Q_inf*sin(alpha) + wake_velocity[0],-Q_inf*cos(alpha) + wake_velocity[1]],n)

		#RHS[4] should be sum of the four foil vortices

		RHS[4] = gamma_foil_list[k]

		gamma_j = solve(A,RHS)
		gamma_j = sum(gamma_j)


		wake_counter += 1 #add another wake vortex
				#move collocation,vortex, LE and TE locations left by Q_inf*dt 


		x_colloc_locs = [i - Q_inf*sin(alpha)*dt for i in x_colloc_locs]
		vortex_x_locs = [i - Q_inf*sin(alpha)*dt for i in vortex_x_locs]
		LE_loc = [LE_loc[0] - Q_inf*dt,LE_loc[1]]
		TE_loc = [TE_loc[0] - Q_inf*dt,TE_loc[1]]
		
		new_wake_loc = [TE_loc[0] + 0.25*Q_inf*dt*cos(alpha),TE_loc[1] - 0.25*Q_inf*dt*sin(alpha)]
		
		vortex_x_locs.append(new_wake_loc[0])
		vortex_z_locs.append(new_wake_loc[1])




	plot(x_colloc_locs,z_colloc_locs,vortex_x_locs,vortex_z_locs,TE_loc,LE_loc) #plots latest location of airfoil, after all iterations 

	# for i in wake_locs:
	# 	plt.plot(i[0],i[1],'r*')

	# plt.ylim((-0.25,0.25))
	# plt.legend()
	# plt.title("Fixed wake vortices and airfoil location after {} iterations".format(ITs))
	# plt.grid()
	# plt.legend()
	plt.show()