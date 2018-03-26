#!/usr/bin/python
import numpy as np 
import scipy


def log_like( sigma_data ):

	sigma_data2 = np.power(sigma_data,2)

	data = np.loadtxt('./data.dat')
	data = data[:,1]

	model = np.loadtxt('./result.txt')
	model = model[:,1]

	N = data.shape[0]

	if N == 1 :
		if np.isnan((model[1])):
			return -1e12
		dif = data[1] - model[1]
		ssn = np.power(dif,2)/sigma_data2

	else:
		if np.isnan(sum(model[:])):
			return -1e12
		dif = data - model
		dif = dif*dif
		ssn = np.sum( dif/sigma_data2 );


	fval =  - 0.5*N*np.log(2.*np.pi) - 0.5*N*np.log(sigma_data2)  - 0.5*ssn

	return fval

#--------------------------------------------------------------------


import sys
import numpy as np

p = np.loadtxt( sys.argv[1] )

N = int(float(sys.argv[2]))


res = log_like( p[N-1] )

np.savetxt('loglike.txt', [res] )
