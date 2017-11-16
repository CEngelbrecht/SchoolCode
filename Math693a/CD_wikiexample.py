import scipy.sparse as sparse
import numpy as np
from numpy import linalg
import numpy as np 
import matplotlib.pyplot as plt

x_0 = np.array([2,1])

b = np.array([1,2])

A = [[4,1],[1,3]]

r = b - A*x_0

print r