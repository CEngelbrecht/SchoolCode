def UMSTOP(n,x_c,x_plus,f,g,typf,retcode,gradtol,steptol,itncount,itnlimit,maxtaken,consecmax):


	termcode = 0 

	if retcode == 1: 
		termcode = 3
	#2b
	elif (max([float(abs(g[i]) * abs(x_plus[i])/max(abs(f),typf)) for i in range(0,n)])) <= gradtol: 
		termcode = 1
	#2c
	elif max([float(abs((x_plus[i] - x_c[i])) / abs(x_plus[i])) for i in range(0,n)]) <= steptol: 

		termcode = 2 

	elif itncount >= itnlimit: 

		termcode = 4

	elif maxtaken == True: 

		consecmax += 1

		if consecmax == 5: 

			termcode = 5
	else: 

		consecmax = 0

	#print("termcode = {} from UMSTOP".format(termcode))

	return termcode,consecmax