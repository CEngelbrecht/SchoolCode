def FDGRAD(n,x_c,f_c,FN,eta):
	
	from numpy import sqrt,zeros,sign,reshape,cbrt
	from copy import copy 

	#sqrteta = sqrt(eta)

	g_c = zeros(n)

	g_c = reshape(g_c,(n,1))


	for j in range(0,n):

		stepsizej = sqrt(eta) * abs(x_c[j]) * sign(x_c[j])

		#stepsizej = cbrt(eta) * abs(x_c[j]) * sign(x_c[j])

		tempj = copy(x_c[j])

		x_c[j] = x_c[j] + stepsizej

		stepsizej = x_c[j] - tempj

		f_j = FN(n,x_c)

		g_c[j] = (f_j - f_c)/stepsizej

		x_c[j] = tempj

	return g_c