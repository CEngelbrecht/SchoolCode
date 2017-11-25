def UMSTOP(n,x_c,x_plus,f,g,typf,retcode,gradtol,steptol,itncount,itnlimit,maxtaken):


	termcode = 0 

	if retcode == 1: 
		termcode = 3

	elif (max([float(abs(g_0[i]) * abs(x_0[i])/max(abs(f),typf)) for i in range(0,n)])) <= gradtol: 
		termcode = 1

	elif max([float((x_plus[i] - x_0[i]) / x_plus[i]) for i in range(0,n)]) <= steptol: 

		termcode = 2 

	elif itncount >= itnlimit: 

		termcode = 4

	elif maxtaken == True: 

		consecmax += 1

		if consecmax == 5: 

			termcode = 5
	else: 

		consecmax = 0

	print("termcode = {}".format(termcode))

	return termcode