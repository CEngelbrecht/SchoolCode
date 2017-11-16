import matplotlib.pyplot as plt 
import numpy as np 
import sympy as sp 

x1,x2 = sp.symbols('x1 x2')

func = x1**2 + x2**2/4 + 4*(x1 - x2)**2 * sp.sin(x2)**2 

f = sp.lambdify((x1,x2),func,"numpy") #A lambdified version: speeds it up, faster than evalf(). Evaluate this with f(x1,x2)

delta = 0.025
x = np.arange(-2.0, 2.0, delta)
y = np.arange(-2.0, 2.0, delta)
X, Y = np.meshgrid(x, y)

Z = X**2 + (Y**2)/4.0 + 4*(X - Y)**2 * np.sin(y)**2

plt.figure()
cp = plt.contour(X, Y, Z)
plt.show()