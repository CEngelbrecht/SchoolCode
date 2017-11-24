'''
Christian Engelbrecht
Math693a Final "Default" Project
'''



#import necessary functions 
import numpy as np 
import sys 

from UMSTOP0 import UMSTOP0
from MODELHESS import MODELHESS
import CHOLDECOMP
import CHOLSOLVE 
from MACHINEPS import MACHINEPS
import UMSTOP
from UMINCK import UMINCK
from FN import FN
from GRAD import GRAD
from HESS import HESS
''' 

Currently testing with rosenbrock function, with n = 2, i.e. f(x) = 100(x2 - x1**2)**2  + (1 - x2)**2
Syntax notes: x_c -> x_k, f_c -> f_k for consistency with previous assignments 
'''

analgrad = True
analhess = True
cheapf = False
factsec = False 
gradtol = 1E-8

n = 2
x_0 = (1.2,1.2) #inital starting tuple
x_0 = np.reshape(np.array(x_0),(len(x_0),1)) #numpy array of (x_0)T


#				Start initialization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
machineps = MACHINEPS()
termcode = UMINCK(n,machineps,x_0)

if termcode < 0: 

	print("The size the of the problem is too small! Exiting now")
	sys.exit()

itncount = 0 

f_0 = FN(n,x_0)

#6
if analgrad == True: 

	g_0 = GRAD(n,x_0)

else: 
	
	pass #call FDGRAD

#7
termcode = UMSTOP0(n,x_0,f_0,g_0,gradtol) #this doesn't work yet 

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

	if not factsec:
		
		L = MODELHESS(n,machineps,H_c)

