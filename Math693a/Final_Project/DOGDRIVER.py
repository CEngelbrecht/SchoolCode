def DOGDRIVER(n,x_c,f_c,FN,g_c,L,s_N,maxstep,steptol,delta): 

	from DOGSTEP import DOGSTEP 
	from TRUSTREGUP import TRUSTREGUP

	from numpy.linalg  import norm 

	retcode = 4

	firstdog = True

	Newtlen = norm(s_N)

	while retcode > 2:
		print "retcode in while loop = {}".format(retcode)

		if firstdog == True:

			try: 
				#accepts the Newton step

				delta, s,Newttaken = DOGSTEP(n,g_c,L,s_N,Newtlen,maxstep,delta,firstdog)
			except ValueError as e: 

				delta,firstdog,Cauchylen,eta,s_SD,v_hat,s,Newttaken = DOGSTEP(n,g_c,L,s_N,Newtlen,maxstep,delta,firstdog)

		elif firstdog == False: 

			try:

				delta, s,Newttaken = DOGSTEP(n,g_c,L,s_N,Newtlen,maxstep,delta,firstdog,Cauchylen = Cauchylen, eta = eta, s_SD =  s_SD,v_hat = v_hat)
			except ValueError as e: 

				delta, firstdog,Cauchylen,eta,s_SD,v_hat,s,Newttaken = DOGSTEP(n,g_c,L,s_N,Newtlen,maxstep,delta, firstdog, Cauchylen = Cauchylen, eta = eta, s_SD =  s_SD,v_hat = v_hat)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Calling Trustregup~~~~~~~~~~~~~~~~~~~

		if retcode == 3:

			try:
				
				delta,retcode,x_plus_prev,f_plus_prev,x_plus,f_plus,maxtaken = TRUSTREGUP(n,x_c,f_c,FN,g_c,L,s,Newttaken,maxstep,steptol,2,L,delta,retcode,x_plus_prev = x_plus_prev,f_plus_prev = f_plus_prev)

			except (ValueError,UnboundLocalError) as e: 

				print e
				delta,retcode,x_plus_prev,f_plus_prev,maxtaken= TRUSTREGUP(n,x_c,f_c,FN,g_c,L,s,Newttaken,maxstep,steptol,2,L,delta,retcode,x_plus_prev = x_plus_prev,f_plus_prev = f_plus_prev)
		else:

			try:

				delta, retcode,x_plus,f_plus,maxtaken = TRUSTREGUP(n,x_c,f_c,FN,g_c,L,s,Newttaken,maxstep,steptol,2,L,delta,retcode)

			except ValueError as e: 


				delta,retcode,x_plus_prev,f_plus_prev,x_plus,f_plus,maxtaken = TRUSTREGUP(n,x_c,f_c,FN,g_c,L,s,Newttaken,maxstep,steptol,2,L,delta,retcode)

	return delta,retcode,x_plus,f_plus,maxtaken