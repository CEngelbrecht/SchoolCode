import numpy as np 

A = np.array([[4,12,-16],[12,37,-43],[-16,-43,98]])

L = np.zeros(A.shape)

n = 3 

for j in range(len(A)):

	L[j][j] = A[j][j] - sum([L[j][i]**2 for i in range(0,j)])
	
	minljj = 0 
	
	for i in range(j+1,n):

		L[i][j] = A[j][i] - sum([L[i][k] * L[j][k] for k in range(0,(j - 1))])
		
		minljj = max(abs(L[i][j]),minljj)
		
	if L[j][j] > minljj**2:
		
		L[j][j] = np.sqrt(L[j][j])
	
	for i in range(j+1,n):
	
		L[i][j] = L[i][j]/L[j][j]
	

print L
		


