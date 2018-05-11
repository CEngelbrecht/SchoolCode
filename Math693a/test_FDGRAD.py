from FDGRAD import FDGRAD 
from Rosenbrock_Driver import GRAD_rosenbrock as GRAD
from Rosenbrock_Driver import FN_rosenbrock as FN
import numpy as np
from MACHINEPS import MACHINEPS

eta = MACHINEPS()

n = 2 

x_0 = np.array([[-1.2],[1.0]])

f_c = FN(n,x_0)

grad = GRAD(n,x_0)

g = FDGRAD(n,x_0,f_c,FN,eta)

print grad
print g

print grad - g 