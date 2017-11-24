def MACHINEPS():
	'''Determines Machine Epsilon'''

	machineps = 1 

	while (1 + machineps) != 1: 

		machineps = machineps/2.0

	machineps = 2 * machineps

	return machineps