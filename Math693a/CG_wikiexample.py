import scipy.sparse as sparse
import numpy as np
from numpy import linalg
import numpy as np 
import matplotlib.pyplot as plt

x = np.array([2,1]) #Inital guess, this is x_0 

b = np.array([1,2])

A = np.array([[4,1],[1,3]])

r = np.dot(A,x) - b 

p = -1*r 

while np.linalg.norm(r) > 1E-8:
	
	#print r,np.linalg.norm(r)

	alpha_num = np.dot(r,r)
	alpha_denom = (np.dot(np.transpose(p),np.dot(A,p)))

	alpha = float(alpha_num)/alpha_denom

	x  = x + alpha * p #Now we have the first step done. Time to find the residual, update p, and move on. 

	r_kplus1  = r + alpha * np.dot(A, p)

	beta = (np.dot(r_kplus1,r_kplus1))/(np.dot(r,r))

	p_kplus1 = -1*r_kplus1 + beta * p

	r = r_kplus1

	p = p_kplus1

print x

