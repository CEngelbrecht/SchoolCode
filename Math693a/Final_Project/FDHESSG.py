def FDHESSG(n,x_c,g_c,GRAD,eta):

	from FDJAC import FDJAC
	from numpy import transpose

	H = FDJAC(n,x_c,g_c,GRAD,eta)

	#print("H as returned by FDJAC = \n{}".format(H))

	

	return (H + transpose(H))/2.0
