def FDJAC(n,x_c,F_c,FVEC,eta):

	from numpy import sqrt,sign,zeros,cbrt
	from copy import copy 

	sqrteta = sqrt(eta)

	J = zeros((n,n))

	stepsizej = sqrt(eta)

	for j in range(0,n):

		stepsizej = sqrt(eta) * abs(x_c[j]) * sign(x_c[j])

		tempj = copy(x_c[j]) #assign current value of x_c[j] to a temp value

		x_c[j] = x_c[j] + stepsizej #perturb x_c + stepsize 

		F_j = FVEC(n,(x_c))

		for i in range(0,n):
			J[i][j] = (F_j[i] - F_c[i])/stepsizej

		x_c[j] = tempj

	return J 
