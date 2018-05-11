from FDHESSF import FDHESSF
from Rosenbrock_Driver import FN_rosenbrock as FN
from Rosenbrock_Driver import HESS_rosenbrock as HESS
import numpy as np 
from MACHINEPS import MACHINEPS

eta = MACHINEPS()

n = 2 

x_0 = np.array([[-1.2],[1.0]])

Hess = HESS(n,x_0)

f_c = FN(n,x_0)

H = FDHESSF(n,x_0,f_c,FN,eta)

print Hess
print H