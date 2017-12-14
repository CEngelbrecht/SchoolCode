def DOGSTEP(n,g_c,L,s_N,Newtlen,maxstep,delta,firstdog,**kwargs): 

	from numpy import dot,transpose,sqrt 
	from numpy.linalg import norm

	#Retrieve kwargs, if any
	if 'Cauchylen' in kwargs: 
		Cauchylen = kwargs["Cauchylen"]
	if 'eta' in kwargs:
		eta = kwargs['eta']
	if 's_SD' in kwargs:
		s_SD = kwargs["s_SD"]
	if 'v_hat' in kwargs:	
		v_hat = kwargs['v_hat']

	if Newtlen <= delta: 

		Newttaken = True

		s = s_N

		delta = Newtlen

		print "Returning the Newton Step"

		return delta,s,Newttaken

	else: 

		Newttaken = False

		if firstdog == True: 

			firstdog = False 

			alpha = norm(g_c)**2

			beta = 0 

			for i in range(0,n): 

				temp = sum([L[j][i] * g_c[j] for j in range(i,n)])

				beta = beta + temp*temp

			s_SD = - (alpha/beta) * g_c

			Cauchylen = (alpha * alpha**(1/2))/beta

			eta = 0.2 + (0.8 * alpha**2 /(beta * abs(dot(transpose(g_c),s_N))))

			v_hat = eta * s_N - s_SD

			if delta == -1: 

				delta = min(Cauchylen,maxstep)

		if (eta * Newtlen) <= delta:

			s = (delta/Newtlen)*s_N

			#return from here
		elif Cauchylen >= delta: #take steepest descent direction

			s = (delta/Cauchylen) * s_SD

			#return from here

		else: 

			temp = dot(transpose(v_hat),s_SD)
			tempv = dot(transpose(v_hat),v_hat)

			Lambda = (-temp + sqrt(temp**2 - tempv * (Cauchylen**2 - delta**2)))*tempv #Check the final tempv term... might not be multiplied

			s = (s_SD + Lambda*v_hat)

		return delta,firstdog,Cauchylen,eta,s_SD,v_hat,s,Newttaken	
