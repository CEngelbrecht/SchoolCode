'''
Christian Engelbrecht
Math693a Final "Default" Project
'''

def UMDRIVER(n,x_0,FN,GRAD,HESS,globalstrat,analgrad,analhess,cheapf,factsec,gradtol,steptol,maxstep,itnlimit,**kwargs):

	import numpy as np 
	import sys 
	from UMSTOP0 import UMSTOP0
	from MODELHESS import MODELHESS
	from CHOLDECOMP import CHOLDECOMP
	from CHOLSOLVE import CHOLSOLVE
	from LINESEARCH import LINESEARCH
	from MACHINEPS import MACHINEPS
	from UMSTOP import UMSTOP
	from UMINCK import UMINCK
	from PLOTTER import PLOTTER
	from DOGDRIVER import DOGDRIVER 
	from FDHESSG import FDHESSG
	from FDGRAD import FDGRAD
	from FDHESSF import FDHESSF

	#print("UMDRIVER called with kwargs {}".format(kwargs))

	#assign kwargs if passed on from calling function
	if 'typf' in kwargs:
		typf = kwargs['typf']
	if 'plot_results' in kwargs:
		plot_results = kwargs['plot_results']
	if 'delta' in kwargs: 
		delta = kwargs['delta']
	if 'fdigits' in kwargs:
		fdigits = kwargs['fdigits']




	#				Book keeping 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	x_list = []
	f_list = []
	g_list = []
	s_list = []
	delta_list = []

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	#				Start initialization
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	machineps = MACHINEPS()
	termcode,eta = UMINCK(n,machineps,x_0,fdigits = fdigits) #needs beefing up
	consecmax = 0.0

	if termcode < 0: 

		x_f = x_0
		print("The size the of the problem is too small! Exiting now")
		sys.exit()

	itncount = 0 

	f_c = FN(n,x_0)

	#6
	if analgrad == True: 
		g_c = GRAD(n,x_0)

	else: 
		g_c = FDGRAD(n,x_0,f_c,FN,eta)

	#7
	termcode = UMSTOP0(n,x_0,f_c,g_c,gradtol) 

	if termcode > 0: 
		print "x_0 is a solution already"
		sys.exit()

	else: 
		if analhess == True: 
			H_c = HESS(n,x_0)

		elif (analgrad and cheapf): 
			H_c =  FDHESSG(n,x_0,g_c,GRAD,eta)

		elif cheapf: 
			H_c = FDHESSF(n,x_0,f_c,FN,eta)

		elif factsec: 
			pass 

		elif (factsec == False): 
			pass 

	x_c = x_0 
	H_0 = H_c
						#Book keeping
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	x_list.append(x_c)
	g_list.append(g_c)
	f_list.append(f_c)

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
						#Done initializing, start iterating 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	while termcode == 0:

		print("UMDRIVER itncount = {}, x_c = \n{} g_c = \n{} H_c = \n{} ".format(itncount,x_c,g_c,H_c))
		itncount += 1 

		if not factsec:
			H_c,L_c = MODELHESS(n,machineps,H_c) #returns model hessian's lower triangular parts(?)

		s_N = CHOLSOLVE(n, g_c, L_c) #s_N is newton search direction  

		if globalstrat == 1:
			#do linesearch 
			retcode,x_plus,f_plus,maxtaken = LINESEARCH(n,x_c,f_c,FN,g_c,s_N,maxstep,steptol) # 0 retcode = found good x_plus 

		elif globalstrat == 2:
			#do hookdriver 

			pass

		elif globalstrat == 3:
			#Call dogdriver 

			delta,retcode,x_plus,f_plus,maxtaken = DOGDRIVER(n,x_c,f_c,FN,g_c,L_c,s_N,maxstep,steptol,delta)

		else: 
			#call linsearch mod 
			pass 

		#10.5
		if globalstrat != 4: 

			if analgrad == True: 

				g_plus = GRAD(n,x_plus)

			else: 
				#call FDGRAD
				g_plus = FDGRAD(n,x_c,f_c,FN,eta)
		#10.6
		termcode,consecmax = UMSTOP(n,x_c,x_plus,f_plus,g_plus,typf,retcode,gradtol,steptol,itncount,itnlimit,maxtaken,consecmax)

		if termcode > 0: 
			#Found a final candidate
			x_f = x_plus

		else:
			#The search continues...
			if analhess == True: 

				H_c = HESS(n,x_plus)

			elif (analgrad == True and cheapf == True):
				#call FDHESSG
				H_c =  FDHESSG(n,x_c,g_c,GRAD,eta)
				

			elif cheapf == True:
				#call FDHESSF
				H_c = FDHESSF(n,x_plus,f_plus,FN,eta)

			elif factsec == True:
				#call BFGSFAC
				pass

			else:
			#call BFGSUNFAC
				pass

			x_c = x_plus
			f_c = f_plus
			g_c = g_plus

			#~~~~~~~~~~~Bookkeeping~~~~~~~~~~~~~~~~
			x_list.append(x_c)
			g_list.append(g_c)
			f_list.append(f_c)
			s_list.append(s_N)
			delta_list.append(delta)

			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	if plot_results:
		PLOTTER(f_list,globalstrat = globalstrat, tag = 'function',)
		PLOTTER(g_list,globalstrat = globalstrat,tag = 'gradient',)
	if plot_results and n == 2 and globalstrat == 1 :
		PLOTTER(x_list,globalstrat = globalstrat,tag = 'direction')
	elif plot_results and n == 2 and globalstrat == 3: 
		PLOTTER(x_list,globalstrat = globalstrat,tag = 'direction',delta_list =  delta_list)

	return x_f,termcode,delta_list
