def FDHESSF(n,x_c,f_c,FN,eta):

	from numpy import cbrt,sign,zeros
	from copy import copy 

	cuberteta = cbrt(eta)

	H = zeros((n,n))

	stepsize = zeros(n)
	fneighbor = zeros(n)

	for i in range(0,n):

		stepsize[i] = cuberteta * abs(x_c[i]) * sign(x_c[i])

		tempi = copy(x_c[i])

		x_c[i] = x_c[i] + stepsize[i]

		fneighbor[i] = FN(n,x_c)

		x_c[i] = tempi

	for i in range(0,n):
		#print i

		tempi = copy(x_c[i])

		x_c[i] = x_c[i] + 2 *stepsize[i]

		fii = FN(n,x_c)

		H[i][i] = ((f_c - fneighbor[i]) + (fii - fneighbor[i]))/(stepsize[i]*stepsize[i])

		x_c[i] = tempi + stepsize[i] 

		for j in range(i,n):
			#print i,j

			tempj = copy(x_c[j])

			x_c[j] = x_c[j] + stepsize[j]

			fij = FN(n,x_c)

			H[i][j] = ((f_c - fneighbor[i]) + (fij - fneighbor[j]))/(stepsize[i]*stepsize[j])

			x_c[j] = tempj

		x_c[i] = tempi


	return H 