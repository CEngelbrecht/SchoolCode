import numpy as np
from UMDRIVER import UMDRIVER
from MACHINEPS import MACHINEPS
import time
import matplotlib.pyplot as plt  

function = "Rosenbrock"

if function == "Rosenbrock":

	from Rosenbrock_Driver import FN_rosenbrock as FN
	from Rosenbrock_Driver import GRAD_rosenbrock as GRAD
	from Rosenbrock_Driver import HESS_rosenbrock as HESS

	n = 2

	x_0 = np.zeros(n)

	x_0 = np.reshape(x_0,(n,1))

	for i in range(0,n,2):

		x_0[i] = -1.2
		x_0[i+1] = 1.0

elif function == "Powell":

	from Powell_Driver import FN_Powell as FN
	from Powell_Driver import GRAD_Powell as GRAD
	from Powell_Driver import HESS_Powell as HESS

	n = 12

	x_0 = np.zeros(n)

	x_0 = np.reshape(x_0,(n,1))

	for i in range(0,n,4):

		x_0[i] = -3.0
		x_0[i+1] = -1.0
		x_0[i+2] = 0 	
		x_0[i+3] = 1.0

globalstrat = 1
analgrad = True
analhess = True
factsec = False
cheapf = True 
fdigits = -1 
typf = 1 
plot_results = True

machineps = MACHINEPS()

gradtol = machineps**(1.0/3.0)
steptol = machineps**(2.0/3.0)
maxstep = 10
itnlimit = 100
delta = 1.0

time_start = time.time()

x_f,termcode,delta_list = UMDRIVER(n,x_0,FN,GRAD,HESS,globalstrat,analgrad,analhess,cheapf, \
						factsec,gradtol,steptol,maxstep,itnlimit,typf = typf,delta = delta,plot_results = plot_results,\
								fdigits = fdigits)
time_end = time.time()

#print("Time taken in {}".format(time_end - time_start))

#print "UMDRIVER terminated with  x_f = {}".format(x_f)

if termcode  == 1: 
	print("termcode = 1. Norm of scaled gradient less than gradtol; x_plus probably is an approximate local minimizer of f(x)")
elif termcode == 2: 
	print("termcode = 2. Scaled distance between last two steps less than steptol, x_plus may be an approximate local minimizer of fx), but it is also possible that the algorithm is making very slow progress and is not near a minimizer. Or,steptol is too large ")
elif termcode == 3: 
	print("termcode = 3. Last global step failed to locate a lower point than x_c. Either x_c is an approximate local minimizer and no more \
accuracy is possible, or an innacurately coded analytic gradient is being used, or the finite difference approximation \
is too innacurate, or steptol is too large")
elif termcode == 4: 
	print("Termcode = 4. Iteration limit exceeded")
elif termcode == 5: 
	print("Termcode = 5. Five consecutive steps of length maxstep have been taken: either f(x) is unbounded below, or f(x) has\
a finite asymptote in some direction, or maxstep is too small")
else: 
	print("Not a recognized termcode")

correct_x = np.zeros(n)

correct_x = np.reshape(correct_x,(n,1))

for i in range(len(correct_x)):
	correct_x[i] = 0

ave = np.average(abs(x_f - correct_x))
print "{:.5E}".format(ave)