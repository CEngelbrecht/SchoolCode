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

machineps = MACHINEPS()

gradtol = machineps**(1.0/3.0)
steptol = machineps**(2.0/3.0)
maxstep = 1E3
itnlimit = 15000
delta = -1.0

UMDRIVER(n,x_0,FN,GRAD,HESS,globalstrat,analgrad,analhess,cheapf,factsec,*typx,*typf,*fdigits,*gradtol,steptol,maxstep,itnlimit,delta)

