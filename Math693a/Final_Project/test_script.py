import numpy as np
from UMDRIVER import UMDRIVER
from MACHINEPS import MACHINEPS
from FN import FN 
from GRAD import GRAD
from HESS import HESS 

'''Watch out for pass by reference https://www.tutorialspoint.com/python/python_functions.htm '''

x_0 = (1.2,1.2) #inital starting tuple
x_0 = np.reshape(np.array(x_0),(len(x_0),1)) #numpy array of (x_0)T
n = len(x_0)


globalstrat = 1 
analgrad = True
analhess = True
factsec = False
cheapf = False 
fdigits = -1 
typf = 1 

machineps = MACHINEPS()

gradtol = machineps**(1.0/3.0)
steptol = machineps**(2.0/3.0)
maxstep = 1E3
itnlimit = 15000
delta = -1.0

x_f,termcode = UMDRIVER(n,x_0,FN,GRAD,HESS,globalstrat,analgrad,analhess,cheapf, \
						factsec,gradtol,steptol,maxstep,itnlimit,typf = typf,delta = delta)

print "UMDRIVER terminated with  x_f = {}".format(x_f)

if termcode  == 1: 
	print("termcode = 1. Norm of scaled gradient less than gradtol; x_plus probably is an approximate local minimizer of f(x)")
elif termcode == 2: 
	print("termcode = 2. Scaled distance between last two steps less than steptol, x_plus may be an approximate local minimizer of fx), but it is also possible that the algorithm is making very slow progress and is not near a minimizer. Or,steptol is too large ")
elif termcode == 3: 
	print("last global step failed to locate a lower point than x_c. Either x_c is an approximate local minimizer and no more \
accuracy is possible, or an innacurately coded analytic gradient is being used, or the finite difference approximation \
is too innacurate, or steptol is too large")
elif termcode == 4: 
	print("Iteration limit exceeded")
elif termcode == 5: 
	print("five consecutive steps of length maxstep have been taken: either f(x) is unbounded below, or f(x) has\
a finite asymptote in some direction, or maxstep is too small")
else: 
	print("Not a recognized termcode")