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

	print("UMDRIVER called with kwargs {}".format(kwargs))

	if 'typf' in kwargs:
		typf = kwargs['typf']

	#				Start initialization
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	machineps = MACHINEPS()
	termcode = UMINCK(n,machineps,x_0) #needs beefing up

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
		pass #call FDGRAD

	#7
	termcode = UMSTOP0(n,x_0,f_c,g_c,gradtol) 

	if termcode > 0: 
		print "x_0 is a solution already"
		sys.exit()

	else: 
		if analgrad == True: 
			H_c = HESS(n,x_0)

		elif (analgrad and cheapf): 
			pass 

		elif cheapf: 
			pass 

		elif factsec: 
			pass 

		elif (factsec == False): 
			pass 

	x_c = x_0 
	H_0 = H_c
						#Done initializing, start iterating 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	while termcode == 0: 

		print("UMDRIVER itncount = {}".format(itncount))
		itncount += 1 

		if not factsec:
			L_c = MODELHESS(n,machineps,H_c) #returns 

		s_N = CHOLSOLVE(n, g_c, L_c)

		if globalstrat == 1: #do linesearch 	
			print("linsearch called with x_c = {}, s_N = {}".format(x_c,s_N))
			retcode,x_plus,f_plus,maxtaken = LINESEARCH(n,x_c,f_c,g_c,s_N,maxstep,steptol) # 0 retcode = found good x_plus 

		elif globalstrat == 2: #do hookdriver 
			pass
		elif globalstrat ==3: #do dogdriver
			pass 
		else: 
			#call linsearch mod 
			pass 

		#10.5
		if globalstrat != 4: 

			if analgrad == True: 

				g_plus = GRAD(n,x_plus)

			else: 
				#call FDGRAD
				pass
		#10.6
		termcode = UMSTOP(n,x_c,x_plus,f_plus,g_plus,typf,retcode,gradtol,steptol,itncount,itnlimit,maxtaken)

		if termcode > 0: 
			#Found a final candidate
			x_f = x_plus

		else:
			#The search continues...
			if analhess == True: 

				H_c = HESS(n,x_plus)

			elif (analgrad == True and cheapf == True):
				#call FDHESSG
				pass 

			elif cheapf == True:
				#call FDHESSF
				pass 

			elif factsec == True:
				#call BFGSFAC
				pass

			else:
			#call BFGSUNFAC
				pass

			x_c = x_plus
			f_c = f_plus
			g_c = g_plus

	return x_f,termcode
