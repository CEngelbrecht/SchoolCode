import numpy as np 
from MACHINEPS import MACHINEPS
from Rosenbrock_Driver import HESS_rosenbrock as HESS
from Rosenbrock_Driver import GRAD_rosenbrock as GRAD
from FDHESSG import FDHESSG

eta = MACHINEPS()
n =2 
x_0 = np.array([[-1.2],[1.0]])

g_c = GRAD(n,x_0)

Hess = HESS(n,x_0)

H = FDHESSG(n,x_0,g_c,GRAD,eta)

print("Real Hessian is \n{}".format(Hess))
print("H as returned from FDHESSG = \n{}".format(H))
