def TRUSTREGUP(n,x_c,f_c,FN, g_c,L,s,Newttaken,maxstep,steptol,steptype,H,delta,retcode,**kwargs): 

	from numpy import dot,transpose 
	from numpy.linalg import norm

	if "x_plus_prev" in kwargs: 
		x_plus_prev = kwargs["x_plus_prev"]
	if "f_plus_prev" in kwargs:
		f_plus_prev = kwargs["f_plus_prev"]

	maxtaken = False 
	alpha = 1E-4
	steplen = norm(s)

	x_plus = x_c + s

	f_plus = FN(n,x_plus)

	delta_f = f_plus - f_c

	initslope = dot(transpose(g_c),s)

	print("initslope = {}".format(initslope))

	if retcode != 3:

		f_plus_prev = 0

	if (retcode == 3) and ((f_plus >= f_plus_prev) or (delta_f > alpha * initslope)): 

		retcode = 0 

		x_plus = x_plus_prev

		f_plus = f_plus_prev

		delta = delta/2.0

		#return from here 

	elif delta_f >= alpha*initslope: 
		#f(x_plus) too large 

		rellength = max([s[i]/x_plus[i] for i in range(0,n)])

		if rellength < steptol: 

			#x_plus  - x_c too small
			retcode = 1 
			x_plus = x_c

			#return from here 
		else: 

			retcode = 2 

			delta_temp = (-initslope * steplen)/(2 * (delta_f - initslope))

			if delta_temp < 0.1*delta: 

				delta = 0.1 * delta

			elif delta_temp > 0.5*delta: 

				delta = 0.5 * delta 

			else: 
				delta = delta_temp

			#return from here 

	else: 
		#f(x_plus) sufficiently small 

		deltaf_pred = initslope

		if steptype == 1:
			#hookstep
			pass 

		else:
			#dogstep
			for i in range(0,n): 

				temp = sum([L[j][i] * s[j] for j in range(i)])

				deltaf_pred = deltaf_pred + (temp*temp/2)

		if retcode != 2 and ((abs(deltaf_pred - delta_f) <= 0.1 * abs(delta_f)) or (delta_f <= initslope)) and Newttaken == False and (delta <= 0.99*maxstep):

			retcode = 3 
			x_plus_prev = x_plus
			f_plus_prev = f_plus 
			delta = min(2*delta,maxstep)

		else: 

			retcode = 0
			if steplen > 0.99 * maxstep: 

				maxtaken = True

			if delta_f >= 0.1 * deltaf_pred: 

				delta = delta/2

			elif delta_f <= 0.75 * deltaf_pred: 

				delta = min(2 * delta,maxstep)

	print("retcode in TRUSTREGUP = {}".format(retcode))

	if retcode == 3: 

		return delta,retcode,x_plus_prev,f_plus_prev,x_plus,f_plus,maxtaken

	else:

		return delta,retcode,x_plus,f_plus,maxtaken