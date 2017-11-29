from numpy.linalg import cholesky
from numpy import array
from GRAD import GRAD
from HESS import HESS 
from MACHINEPS import MACHINEPS
from MODELHESS import MODELHESS
from CHOLSOLVE import CHOLSOLVE

machineps = MACHINEPS()

n = 4
x_0 = array([[1.2],[1.2],[1.2],[1.2]])

g_c = GRAD(n,x_0)

H_c = HESS(n,x_0)

L_c = MODELHESS(n,machineps,H_c)

s_N = CHOLSOLVE(n, g_c, L_c)

print s_N
#cholesky_for_realsky = cholesky(H_c)
#print H_c
#print L 
#print cholesky_for_realsky
