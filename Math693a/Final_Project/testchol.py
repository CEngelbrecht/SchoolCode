from numpy.linalg import cholesky
from numpy import array
from HESS import HESS 
from MACHINEPS import MACHINEPS
from MODELHESS import MODELHESS

machineps = MACHINEPS()

n = 4
x_0 = array([[1.2],[1.2],[1.2],[1.2]])

H_c = HESS(n,x_0)

L = MODELHESS(n,machineps,H_c)

H_c = HESS(n,x_0)

cholesky_for_realsky = cholesky(H_c)
print H_c
print L 
print cholesky_for_realsky
