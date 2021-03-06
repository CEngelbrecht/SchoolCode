def LINESEARCH(n,x_c,f_c,FN,g_c,p, maxstep,steptol): 

	from rosenbrock_2Nd_translated import rosenbrock
	from numpy import dot,sqrt,transpose,array

	maxtaken = False

	retcode = 2

	alpha = 1E-4

	Newtlen = 1 

	if Newtlen > maxstep: 

		p = p * (maxstep/Newtlen)

		Newtlen = maxstep

	initslope = dot(transpose(g_c),p) #gTranspose * p *

	rellength = max([(abs(p[i]))/(x_c[i]) for i in range(0,n)])

	minlambda = steptol/rellength

	Lambda = 1.0

	while retcode >= 2: 

		x_plus = x_c + Lambda * p
		f_plus = FN(n,x_plus) #evaluates function at x_plus 

		if f_plus <= f_c + alpha * Lambda * initslope:  #satisfactory x_plus 
			retcode = 0

			if (Lambda == 1.0) and (Newtlen > 0.99 * maxstep): 
				maxtaken == True

		elif Lambda < minlambda: #no satisfactory x_plus can be found sufficiently distinct from x_c
			retcode = 1 
			x_plus = x_c

		else: 

			if Lambda == 1.0: 
				lambdaTemp = -initslope/(2 * (f_plus - f_c - initslope))

			else: 

				Right_Vector = array([[float(f_plus  - f_c - Lambda*initslope)],[float(f_plusprev - f_c - lambdaPrev*initslope)]])
				Left_Matrix = array([[(1.0/(Lambda**2)),-1.0/(lambdaPrev**2)],[-lambdaPrev/(Lambda**2), Lambda/(lambdaPrev**2)]])
				a,b = (1.0/(Lambda - lambdaPrev)) * dot(Left_Matrix,Right_Vector)
				a,b = float(a),float(b) #extracts floats from the returned 1x1 arrays 

				disc = float(b**2 - 3* a * initslope)

				if a == 0: #cubic is a quadratic
					lambdaTemp = -1.0 * initslope/(2 * b)

				else: 
					lambdaTemp = (-b + sqrt(disc))/(3.0 * a)

				if lambdaTemp > 0.5*Lambda:
					lambdaTemp = 0.5*Lambda

			lambdaPrev = Lambda
			f_plusprev = f_plus

			if lambdaTemp <= 0.1*Lambda: 
				Lambda = 0.1*Lambda
			else:
				Lambda = lambdaTemp

	return retcode,x_plus,f_plus,maxtaken
