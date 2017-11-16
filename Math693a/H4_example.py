import scipy.sparse as sparse
import numpy as np
from numpy import linalg
import numpy as np 
import matplotlib.pyplot as plt

A = np.array([1,4,4,9,9,9,16,16,16,16,25,25,25,25,25])

x = np.zeros(len(A))

b = np.ones(len(A))

r = np.dot(A,x) - b 

p = -1 * r 

# while np.linalg.norm(r) > 1E-8:

for k in range(3):

	print r,np.linalg.norm(r)

	alpha_num = np.dot(r,r)
	alpha_denom = (np.dot(np.transpose(p),np.dot(A,p)))

	alpha = float(alpha_num)/alpha_denom

	x  = x + alpha * p #Now we have the first step done. Time to find the residual, update p, and move on. 

	r_kplus1  = r + alpha * np.dot(A, p)

	beta = (np.dot(r_kplus1,r_kplus1))/(np.dot(r,r))

	p = r_kplus1 - beta * p

	r = r_kplus1

print x

